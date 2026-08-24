"""Merges scored calls across callers into one best-call-per-cluster VCF.

The clustering and per-cluster winner-selection logic below is ported
verbatim (in semantics) from the paper's evaluation code (`extract_overlap_ids`,
`reciprocal_overlap`, `compute_reciprocal_overlap`, `compute_sv_clusters`,
`keep_max_qual_dicast`). Two deliberate deviations from the reference:

* `compute_reciprocal_overlap` is called with `df2` being an explicit
  `.copy()` of `df1` rather than the literal same object -- the installed
  bioframe (0.8.0) raises ``ValueError`` when `bf.closest` is given the same
  dataframe object twice (`dicast/multi.py` hit the same issue and worked
  around it with a single-dataframe self-join call; here we keep the
  reference implementation's two-dataframe call shape, just with a copy, so
  the semantics match the eval code exactly).
* `keep_max_qual_dicast` zeroed out every non-winning row's `dicast_qual` in
  place, because the eval needed every row to survive for PR-curve
  accounting. Here we only need the winner, so `_cluster_winner_index`
  returns the winning row's index and the caller keeps only that row.
"""

import os
import re
import logging

import pandas as pd
import bioframe as bf
import networkx as nx


# Thresholds specified by the merge feature: DEL/DUP/INV cluster by
# reciprocal overlap, INS by breakpoint distance.
RECIPROCAL_OVERLAP_THR = 0.5
INS_BREAKPOINT_DISTANCE_THR = 200

# Callers excluded from the "non-population" pool in the winner-selection
# rule below. Kept as the literal two-entry list from the reference
# implementation for fidelity, even though dicast itself only ever uses 'pav'.
POPULATION_CALLERS = ['pav', 'gnomad']

DICAST_THRESHOLD_DEFAULT = 0.4

SYMBOLIC_ALT = {'DEL': '<DEL>', 'DUP': '<DUP>', 'INS': '<INS>', 'INV': '<INV>'}


def extract_overlap_ids(df1: pd.DataFrame, df2: pd.DataFrame, sv_type: str, max_dist_overlap: int = 500,
                        min_len_overlap: float = 0.7, dup_is_ins: bool = False, ignore_len: bool = False,
                        ignore_type: bool = False, single_bp: bool = False) -> pd.DataFrame:
    """ Extract IDs of overlapping SVs from two dataframes (breakpoint-distance matcher).

    Ported verbatim from the paper's evaluation code. Used here for INS
    (breakpoint distance, `ignore_len=True`); DEL/DUP/INV use
    `compute_reciprocal_overlap` instead.

    Args:
        df1 (pd.Dataframe): Pandas dataframe with SVs
        df2 (pd.Dataframe): Pandas dataframe with SVs
        sv_type (str): Type of SVs to consider
        max_dist_overlap (int, optional): Maximum distance between two start and two end breakpoints for an SV to be considered the same. Defaults to 500.
        min_len_overlap (float, optional): Minimum SV length similarity. Defaults to 0.7.

    Returns:
        pd.Dataframe: Dataframe with IDs of overlapping SVs
    """

    # Filter dataframes by sv_type
    if not ignore_type:
        if dup_is_ins and (sv_type == 'INS' or sv_type == 'DUP'):
            df1 = df1[(df1['sv_type'] == 'INS') | (df1['sv_type'] == 'DUP')].copy().reset_index(drop=True)
            df2 = df2[(df2['sv_type'] == 'INS') | (df2['sv_type'] == 'DUP')].copy().reset_index(drop=True)

            # Create a mask to fetch INS with given length
            mask_ins1 = (df1['sv_type'] == 'INS') & (df1['sv_len'].notna())
            df1.loc[mask_ins1, 'end'] = df1.loc[mask_ins1, 'start'] + df1.loc[mask_ins1, 'sv_len']
            mask_ins2 = (df2['sv_type'] == 'INS') & (df2['sv_len'].notna())
            df2.loc[mask_ins2, 'end'] = df2.loc[mask_ins2, 'start'] + df2.loc[mask_ins2, 'sv_len']

        else:
            df1 = df1[df1['sv_type'] == sv_type].copy().reset_index(drop=True)
            df2 = df2[df2['sv_type'] == sv_type].copy().reset_index(drop=True)

    # Remove variants that are artefacts and bring difficulties to the clustering
    df1 = df1[~((df1['start'] == 1) & (df1['sv_len'] > 50000000))].copy().reset_index(drop=True)
    df2 = df2[~((df2['start'] == 1) & (df2['sv_len'] > 50000000))].copy().reset_index(drop=True)

    # Make sure data types are correct
    df1['start'] = df1['start'].astype(int)
    df1['end'] = df1['end'].astype(int)
    df1['sv_len'] = df1['sv_len'].astype(float)
    df2['start'] = df2['start'].astype(int)
    df2['end'] = df2['end'].astype(int)
    df2['sv_len'] = df2['sv_len'].astype(float)

    # Overlap SVs
    finished_overlap = []
    k = 10
    while True:

        # Get closest variants
        closest_intervals = bf.closest(df1, df2, suffixes=('_1', '_2'), k=k)

        # Group by id_1 and get the closest distance
        closest_intervals_grouped = closest_intervals.groupby('id_1').agg({'distance': max})

        # For all variants that have a closest variant with distance > min_dist_overlap, we found all possible overlaps
        finished_ids = closest_intervals_grouped[closest_intervals_grouped['distance'] > max_dist_overlap].index
        finished_overlap.append(closest_intervals[closest_intervals['id_1'].isin(finished_ids)])

        # For all variants that have a closest variant with distance <= min_dist_overlap, we need to find more overlaps
        not_finished_ids = closest_intervals_grouped[closest_intervals_grouped['distance'] <= max_dist_overlap].index

        # If there are no more overlaps to find, we are done
        if len(not_finished_ids) == 0 or k > 100000:
            break

        # Continue with the not finished ids
        df1 = df1[~df1['id'].isin(finished_ids)]
        k *= 2

    closest_intervals = pd.concat(finished_overlap, ignore_index=True)

    closest_intervals = closest_intervals.dropna(subset=['id_1', 'id_2']).reset_index(drop=True)
    closest_intervals['diff_start'] = abs(closest_intervals['start_1'] - closest_intervals['start_2'])
    closest_intervals['diff_end'] = abs(closest_intervals['end_1'] - closest_intervals['end_2'])
    closest_intervals['diff_len'] = closest_intervals.apply(lambda x: min([x['sv_len_1'], x['sv_len_2']]) / max([x['sv_len_1'], x['sv_len_2']]), axis=1)

    # Create a mask to get all entries where sv_len_1 or sv_len_2 is nan
    mask_len_na = (closest_intervals['sv_len_1'].isna() | closest_intervals['sv_len_2'].isna())

    # Create a mask to get all entries with diff_start < max_dist_overlap and diff_end < max_dist_overlap
    mask_dist_both = (closest_intervals['diff_start'] < max_dist_overlap) & (closest_intervals['diff_end'] < max_dist_overlap)
    mask_dist_single = (closest_intervals['diff_start'] < max_dist_overlap) | (closest_intervals['diff_end'] < max_dist_overlap)
    if single_bp:
        mask_dist = mask_dist_single
    else:
        mask_dist = mask_dist_both

    # Create a mask to get all entries with diff_len > min_len_overlap
    mask_len = (closest_intervals['diff_len'] > min_len_overlap)

    # Apply masks
    if ignore_len:
        overlapping_svs = closest_intervals[mask_dist].copy().reset_index(drop=True)
    else:
        overlapping_svs = closest_intervals[mask_dist & (mask_len | mask_len_na)].copy().reset_index(drop=True)

    # Remove self overlaps
    overlapping_svs = overlapping_svs[overlapping_svs['id_1'] != overlapping_svs['id_2']].copy().reset_index(drop=True)

    # Change back INS end position to original value
    if dup_is_ins:
        overlapping_svs.loc[overlapping_svs['sv_type_1'] == 'INS', 'end_1'] = overlapping_svs.loc[overlapping_svs['sv_type_1'] == 'INS', 'start_1'] + 1
        overlapping_svs.loc[overlapping_svs['sv_type_2'] == 'INS', 'end_2'] = overlapping_svs.loc[overlapping_svs['sv_type_2'] == 'INS', 'start_2'] + 1

    return overlapping_svs.reset_index(drop=True)


def reciprocal_overlap(row: pd.Series) -> float:
    """ Computes reciprocal overlap between two intervals for a row produced by bf.closest.

    Ported verbatim from the paper's evaluation code.

    Args:
        row (pd.Series): Row in a dataframe

    Returns:
        float: Reciprocal overlap
    """

    A_start, A_end, B_start, B_end = row['start_1'], row['end_1'], row['start_2'], row['end_2']
    overlap_size = max(0, min(A_end, B_end) - max(A_start, B_start))
    size_A = A_end - A_start
    size_B = B_end - B_start

    if size_A == 0 or size_B == 0 or overlap_size <= 0:
        return 0

    overlap_A = overlap_size / size_A
    overlap_B = overlap_size / size_B

    return min(overlap_A, overlap_B)


def compute_reciprocal_overlap(df1: pd.DataFrame, df2: pd.DataFrame, sv_type: str, overlap_threshold: float) -> pd.DataFrame:
    """ Computes reciprocal overlap between two dataframes (DEL/DUP/INV matcher).

    Ported verbatim from the paper's evaluation code. `df2` must not be the
    same object as `df1` -- pass a `.copy()` (see module docstring).

    Args:
        df1 (pd.Dataframe): Pandas dataframe with SVs
        df2 (pd.Dataframe): Pandas dataframe with SVs
        sv_type (str): Type of SVs to consider
        overlap_threshold (float): Minimum reciprocal overlap.

    Returns:
        pd.Dataframe: Dataframe with IDs of overlapping SVs
    """

    # Filter dataframes by sv_type
    df1 = df1[df1['sv_type'] == sv_type].copy().reset_index(drop=True)
    df2 = df2[df2['sv_type'] == sv_type].copy().reset_index(drop=True)

    closest_intervals = bf.closest(df1, df2, suffixes=('_1', '_2'), k=10)  # k=10 to get all overlapping variants
    closest_intervals = closest_intervals.dropna(subset=['id_1', 'id_2']).reset_index(drop=True)

    closest_intervals['reciprocal_overlap'] = closest_intervals.apply(reciprocal_overlap, axis=1)
    closest_intervals = closest_intervals[closest_intervals['reciprocal_overlap'] > overlap_threshold].copy().reset_index(drop=True)
    closest_intervals = closest_intervals[closest_intervals['id_1'] != closest_intervals['id_2']].copy().reset_index(drop=True)

    return closest_intervals


def compute_sv_clusters(df: pd.DataFrame, df_overlap: pd.DataFrame, id_column: str) -> pd.DataFrame:
    """ Computes clusters of SVs based on an overlap dataframe (connected components).

    Ported verbatim from the paper's evaluation code.

    Args:
        df (pd.DataFrame): Dataframe with SVs
        df_overlap (pd.DataFrame): Dataframe with overlapping SVs
        id_column (str): Name of the column with the SV IDs

    Returns:
        pd.DataFrame: Dataframe with SVs and cluster labels
    """

    # Create a graph from the overlap dataframe
    df = df.copy()
    G = nx.from_pandas_edgelist(df_overlap, 'id_1', 'id_2')

    # Find connected components. This will return a list of sets.
    components = list(nx.connected_components(G))

    # Create a mapping from node to component
    node_to_component = {}
    for i, component in enumerate(components):
        for node in component:
            node_to_component[node] = i

    # Map the components back to the first dataframe
    df['cluster'] = df[id_column].map(node_to_component)

    # For those without a cluster, assign a new unique cluster label
    next_cluster = len(components)
    for idx, row in df.iterrows():
        if pd.isna(row['cluster']):
            df.at[idx, 'cluster'] = next_cluster
            next_cluster += 1

    df['cluster'] = df['cluster'].astype(int)

    return df


def _cluster_winner_index(group: pd.DataFrame, dicast_threshold: float = DICAST_THRESHOLD_DEFAULT):
    """ Selects the winning row index within one SV cluster.

    Population-aware selection rule ported from the paper evaluation's
    `keep_max_qual_dicast`: among the cluster's rows from non-population
    callers (`caller` not in POPULATION_CALLERS), pick the highest
    `dicast_qual` if any clear `dicast_threshold`; otherwise fall back to
    the highest `dicast_qual` in the whole cluster. Unlike the eval version
    (which zeroed every other row's `dicast_qual` for PR-curve accounting),
    this only returns the winning index -- the caller keeps just that row.

    Args:
        group (pd.DataFrame): One cluster's rows (must retain the index of
            the dataframe the caller will `.loc[]` the winner out of).
        dicast_threshold (float): Population-aware qual threshold.

    Returns:
        Index label of the winning row. If every row's dicast_qual is NaN
        (nothing to rank on), the first row's index is returned.
    """

    non_pop_calls = group[~group['caller'].isin(POPULATION_CALLERS)]
    high_qual_calls = non_pop_calls[non_pop_calls['dicast_qual'] >= dicast_threshold]

    if len(high_qual_calls) > 0:
        return high_qual_calls['dicast_qual'].idxmax()

    if group['dicast_qual'].notna().any():
        return group['dicast_qual'].idxmax()

    # All dicast_qual values in the cluster are NaN -- idxmax() would raise
    # (or, depending on pandas version, silently return NaN); keep the
    # first row deterministically rather than crash or pick a bogus index.
    return group.index[0]


def _clusters_for_sv_type(df: pd.DataFrame, sv_type: str) -> pd.DataFrame:
    """ Clusters one SV type's rows of `df` using the matcher appropriate for it.

    DEL/DUP/INV: reciprocal overlap > 0.5 (compute_reciprocal_overlap).
    INS: breakpoint distance < 200bp, SV length ignored (extract_overlap_ids).

    Args:
        df (pd.DataFrame): All scored calls (any sv_type).
        sv_type (str): SV type to cluster.

    Returns:
        pd.DataFrame: `df`'s rows for this sv_type with a 'cluster' column
            (cluster ids local to this sv_type, 0-based).
    """

    df_type = df[df['sv_type'] == sv_type].copy().reset_index(drop=True)
    if df_type.empty:
        return df_type

    if sv_type == 'INS':
        pairs = extract_overlap_ids(df, df.copy(), sv_type,
                                    max_dist_overlap=INS_BREAKPOINT_DISTANCE_THR, ignore_len=True)
    else:
        pairs = compute_reciprocal_overlap(df, df.copy(), sv_type, RECIPROCAL_OVERLAP_THR)

    return compute_sv_clusters(df_type, pairs, 'id')


def cluster_calls(df: pd.DataFrame) -> pd.DataFrame:
    """ Clusters all of `df`'s rows, per sv_type, into a single frame with
    globally unique 'cluster' ids.

    Args:
        df (pd.DataFrame): Scored calls (id, sv_type, chrom, start, end,
            sv_len, caller, dicast_qual, ... columns).

    Returns:
        pd.DataFrame: `df`'s rows (any sv_type not DEL/DUP/INS/INV is
            dropped -- dicast doesn't score anything else) with a 'cluster'
            column, globally unique across sv_types.
    """

    clustered_parts = []
    next_cluster_offset = 0
    for sv_type in sorted(df['sv_type'].unique()):
        if sv_type not in SYMBOLIC_ALT:
            continue
        part = _clusters_for_sv_type(df, sv_type)
        if part.empty:
            continue
        part = part.copy()
        part['cluster'] = part['cluster'] + next_cluster_offset
        next_cluster_offset = int(part['cluster'].max()) + 1
        clustered_parts.append(part)

    if not clustered_parts:
        empty = df.iloc[0:0].copy()
        empty['cluster'] = pd.Series(dtype=int)
        return empty

    return pd.concat(clustered_parts, ignore_index=True)


def select_merged_calls(df: pd.DataFrame, dicast_threshold: float = DICAST_THRESHOLD_DEFAULT) -> pd.DataFrame:
    """ Clusters all scored calls across callers per SV type and keeps only
    the winning call of each cluster.

    Args:
        df (pd.DataFrame): Scored calls, as read from *.SVs.dicast.tsv.
        dicast_threshold (float): Population-aware qual threshold (see
            `_cluster_winner_index`).

    Returns:
        pd.DataFrame: One row per cluster (same columns as `df`, minus the
            internal 'cluster' column), sorted by chrom then start.
    """

    if df.empty:
        return df.copy()

    clustered = cluster_calls(df)
    if clustered.empty:
        return clustered.drop(columns=['cluster'])

    winner_indices = [
        _cluster_winner_index(group, dicast_threshold)
        for _, group in clustered.groupby('cluster', sort=False)
    ]
    winners = clustered.loc[winner_indices].drop(columns=['cluster'])
    winners = winners.sort_values(['chrom', 'start']).reset_index(drop=True)
    return winners


def genotype_to_gt(genotype) -> str:
    """ Converts a scored TSV's genotype value (e.g. the string "(1, 1)" or
    "(None, None)") into a VCF GT string ("1/1", "./.").

    Args:
        genotype: The 'genotype' column value for one row.

    Returns:
        str: A VCF-style GT string; "./." for anything missing or unparseable.
    """

    if pd.isna(genotype):
        return './.'

    text = str(genotype).strip()
    match = re.match(r'^\(\s*([^,]+?)\s*,\s*([^,]+?)\s*\)$', text)
    if not match:
        return './.'

    alleles = []
    for token in match.groups():
        token = token.strip()
        if token in ('None', 'nan', 'NA', ''):
            alleles.append('.')
        else:
            try:
                alleles.append(str(int(float(token))))
            except ValueError:
                alleles.append('.')
    return '/'.join(alleles)


def read_fai_contigs(fai_path: str) -> list:
    """ Reads an .fai index into a list of (contig_name, length) in file order. """

    contigs = []
    with open(fai_path) as f:
        for line in f:
            if not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            contigs.append((fields[0], int(fields[1])))
    return contigs


def _vcf_header_lines(sample: str, contigs: list) -> list:
    """ Builds the header lines of the minimal merged VCF (fresh header --
    this does not attempt to merge the input VCFs' headers). """

    lines = [
        '##fileformat=VCFv4.2',
        '##source=dicast',
    ]
    for name, length in contigs:
        lines.append(f'##contig=<ID={name},length={length}>')
    lines += [
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">',
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the structural variant">',
        '##INFO=<ID=CALLER,Number=1,Type=String,Description="Caller that produced the winning call for this cluster">',
        '##INFO=<ID=DQ,Number=1,Type=String,Description="Dicast Quality Score">',
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + sample,
    ]
    return lines


def build_merged_vcf(scores_path: str, out_path: str, sample: str, fai_path: str,
                     dicast_threshold: float = DICAST_THRESHOLD_DEFAULT) -> tuple:
    """ Clusters the scored calls in `scores_path` and writes one record per
    cluster (the winner) to a minimal, standalone VCF at `out_path`.

    A fresh, minimal header is written (fileformat/source/contig/INFO/FILTER/
    FORMAT lines) rather than attempting to merge the input VCFs' headers --
    plain text is used instead of vcfpy since the header this needs (a
    handful of INFO/FORMAT lines plus contigs from the .fai) is simpler and
    more robust to hand-write than to coax out of vcfpy's header API for a
    VCF that isn't derived from any single input VCF.

    Args:
        scores_path (str): Path to the sample's *.SVs.dicast.tsv
        out_path (str): Path to write the merged VCF to
        sample (str): Sample name (used as the VCF's single sample column)
        fai_path (str): Reference .fai (defines contig order and lengths)
        dicast_threshold (float): Population-aware qual threshold

    Returns:
        (int, int): (number of merged output records, number of input scored calls)
    """

    df = pd.read_csv(scores_path, sep='\t', low_memory=False,
                     dtype={'sample': str, 'cohort_samples': str})
    input_count = len(df)

    winners = select_merged_calls(df, dicast_threshold=dicast_threshold)

    contigs = read_fai_contigs(fai_path)
    contig_order = {name: i for i, (name, _) in enumerate(contigs)}

    # Calls on contigs the fai doesn't know about shouldn't happen in
    # practice (VariantPrep.filter_variants already restricts to canonical
    # chroms derived from the fai) but are dropped defensively rather than
    # crashing the VCF writer.
    known_mask = winners['chrom'].isin(contig_order)
    dropped = len(winners) - int(known_mask.sum())
    if dropped:
        logging.warning(f'{dropped} winning call(s) on contigs absent from --fai were dropped from the merged VCF')
    known = winners[known_mask].copy()

    known['_contig_rank'] = known['chrom'].map(contig_order)
    known = known.sort_values(['_contig_rank', 'start']).reset_index(drop=True)

    lines = _vcf_header_lines(sample, contigs)
    for _, row in known.iterrows():
        sv_type = row['sv_type']
        alt = SYMBOLIC_ALT.get(sv_type, f'<{sv_type}>')
        end = int(row['end'])
        sv_len = row['sv_len']
        filt = row['filter'] if pd.notna(row['filter']) and str(row['filter']).strip() else '.'
        # Multi-value FILTERs arrive comma-joined from the input parser; the
        # VCF spec separates them with ';'.
        filt = str(filt).replace(', ', ';').replace(',', ';')
        dq = row['dicast_qual']
        dq_str = 'NA' if pd.isna(dq) else str(dq)

        info_parts = [f'SVTYPE={sv_type}', f'END={end}']
        if pd.notna(sv_len):
            info_parts.append(f'SVLEN={int(round(sv_len))}')
        info_parts.append(f'CALLER={row["caller"]}')
        info_parts.append(f'DQ={dq_str}')

        gt = genotype_to_gt(row['genotype'] if 'genotype' in row else None)

        lines.append('\t'.join([
            str(row['chrom']), str(int(row['start'])), str(row['id']), 'N', alt,
            '.', str(filt), ';'.join(info_parts), 'GT', gt,
        ]))

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with open(out_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')

    return len(known), input_count
