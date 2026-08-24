import pandas as pd
import bioframe as bf
import networkx as nx


# Same thresholds used by the old cohort-mode variant matching.
RECIPROCAL_OVERLAP_THR = 0.5
INS_BREAKPOINT_DISTANCE_THR = 200


def _reciprocal_overlap(row: pd.Series) -> float:
    """ Computes reciprocal overlap between two intervals for a row produced by bf.closest. """

    a_start, a_end, b_start, b_end = row['start_1'], row['end_1'], row['start_2'], row['end_2']
    overlap_size = max(0, min(a_end, b_end) - max(a_start, b_start))
    size_a = a_end - a_start
    size_b = b_end - b_start

    if size_a == 0 or size_b == 0 or overlap_size <= 0:
        return 0

    return min(overlap_size / size_a, overlap_size / size_b)


def _matching_pairs(df: pd.DataFrame, sv_type: str) -> pd.DataFrame:
    """ Finds pairs of variants of one SV type that describe the same event.

    DEL/DUP are matched by reciprocal overlap, INS by breakpoint distance
    (mirrors the matching logic the old cohort mode used).
    """

    df_type = df[df['sv_type'] == sv_type].copy().reset_index(drop=True)
    if len(df_type) < 2:
        return pd.DataFrame()

    df_type['start'] = df_type['start'].astype(int)
    df_type['end'] = df_type['end'].astype(int)

    closest = bf.closest(df_type, k=10, suffixes=('_1', '_2'))
    closest = closest.dropna(subset=['gid_1', 'gid_2']).reset_index(drop=True)
    closest = closest[closest['gid_1'] != closest['gid_2']].copy().reset_index(drop=True)

    if sv_type == 'INS':
        is_match = (closest['start_1'] - closest['start_2']).abs() < INS_BREAKPOINT_DISTANCE_THR
    else:
        is_match = closest.apply(_reciprocal_overlap, axis=1) > RECIPROCAL_OVERLAP_THR

    return closest[is_match]


def find_rescue_candidates(own_variant_dfs: dict) -> dict:
    """ Determines, for each sample in a multi-sample group, which variants found by the
    OTHER samples' callers should be rescued (i.e. also scored against this sample's own BAM).

    A rescue candidate for sample S is a variant clustered with a variant from another
    sample and for which S itself has no matching variant among its own caller calls.

    Args:
        own_variant_dfs (dict): Maps sample name to that sample's own (filtered) variant
            dataframe, as produced by VariantPrep for a single sample.

    Returns:
        dict: Maps sample name to a dataframe of rescue candidate rows to add for that
            sample (same columns as the input dataframes, 'sample' and 'caller' rewritten
            to reflect the target sample and the rescue origin).
    """

    samples = list(own_variant_dfs.keys())
    # Empty per-sample frames are dropped before concat (pandas deprecates
    # concatenating empty entries); the fallback keeps the columns when every
    # sample's frame is empty.
    frames = [df for df in own_variant_dfs.values() if not df.empty] or list(own_variant_dfs.values())
    combined = pd.concat(frames, ignore_index=True)
    combined['gid'] = combined.index.astype(str)

    pairs = pd.concat(
        [_matching_pairs(combined, sv_type) for sv_type in combined['sv_type'].unique()],
        ignore_index=True
    )

    graph = nx.Graph()
    graph.add_nodes_from(combined['gid'])
    if not pairs.empty:
        graph.add_edges_from(zip(pairs['gid_1'], pairs['gid_2']))
    clusters = nx.connected_components(graph)

    combined_by_gid = combined.set_index('gid')
    rescue_rows = {sample: [] for sample in samples}

    for cluster in clusters:
        members = combined_by_gid.loc[list(cluster)]
        samples_in_cluster = set(members['sample'])
        missing_samples = [sample for sample in samples if sample not in samples_in_cluster]

        if not missing_samples:
            continue

        # Representative variant to transplant: highest caller-reported qual in the cluster.
        representative = members.sort_values('qual', ascending=False, na_position='last').iloc[0]

        for sample in missing_samples:
            row = representative.copy()
            row['sample'] = sample
            row['id'] = f"rescue_{representative['sample']}_{representative['caller']}_{representative['id']}"
            row['caller'] = f"rescue:{representative['sample']}:{representative['caller']}"
            rescue_rows[sample].append(row)

    return {
        sample: (pd.DataFrame(rescue_rows[sample]).reset_index(drop=True)
                 if rescue_rows[sample] else own_variant_dfs[sample].iloc[0:0].copy())
        for sample in samples
    }
