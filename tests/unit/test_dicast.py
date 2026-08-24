"""Unit tests for the ``dicast.cli`` orchestration module.

Importing it pulls the full dependency closure (vcfpy, pysam, xgboost,
bioframe, networkx, pyBigWig), which is only available in the project env.

Everything heavy is mocked via ``monkeypatch``:

* ``ReferenceAnnotator`` / ``AlignmentAnnotatorIllumina`` / ``Dicast`` (real
  annotation, model and BAM I/O) are replaced by record-the-call stubs on the
  loaded ``dicast`` module. We assert the wrappers call the expected stub
  methods and write the expected outputs -- no real models or BAMs are needed
  for those.

Synthetic inputs (feature TSVs, scores TSV, a vcfpy-parseable VCF) are built in
``tmp_path`` so the tests are fully self-contained. EXPECTED values are always
derived from what the test writes, never from running the function first.

One integration test at the bottom runs the real ``call`` pipeline end-to-end
(no mocks) against the shipped ``tests/data/`` chr21 demo dataset.

Dropped vs. the older lucid/dicast dev line's test_dicast.py
--------------------------------------------------------------
* ``score_variants``'s old BND dispatch case (asserting a
  ``score_translocations`` heuristic call for ``sv_type == 'BND'``) was
  dropped. In this repo, ``score_variants`` only special-cases ``'INV'``;
  every other sv_type -- including 'BND', if ever passed -- goes through the
  XGBoost model-loading branch (``dicast_BND.json``). ``Dicast.score_
  translocations`` still exists on the model class (see
  tests/unit/test_model.py) but ``dicast.py`` no longer calls it from
  anywhere, and no ``BND`` model ships in ``models/`` nor is BND in the
  default ``sv_types``, so that code path is dead from the CLI's
  perspective. See findings below.
* No cohort VCF/CSV mode tests existed in the old module to begin with
  (that logic lived in ``dicast.utils.sample_vcf_to_dataframe``, which
  this repo removed entirely -- confirmed absent from
  ``dicast/utils.py``), so nothing needed dropping there.

New coverage added for this repo
---------------------------------
* ``score_variants``'s ``--pop`` fallback: prefers ``dicast_{sv}_pop.json``
  for DEL/INS when it exists, falls back to the standard model when it does
  not, and never applies to other sv types (e.g. DUP) even with ``--pop``
  set.
* ``combine_feature_files``'s explicit ``dtype={'sample': str,
  'cohort_samples': str}`` on the raw TSV read -- without it, an all-digit
  sample name or a numeric-looking value would be silently re-typed by
  pandas' inference.
* An end-to-end integration test running ``python3 -m dicast call`` on the
  real ``tests/data/`` demo dataset (chr21, 20 real HG002 DEL/INS calls,
  caller label ``delly``).
"""
from __future__ import annotations

import subprocess
import sys
from types import SimpleNamespace

import pandas as pd
import pytest
import vcfpy

from tests.conftest import REPO_DIR
from dicast import cli as dicast


# ---------------------------------------------------------------------------
# Shared helpers.
# ---------------------------------------------------------------------------

class _MethodRecorder:
    """Records every call to a named method (positional + keyword args)."""

    def __init__(self):
        self.calls = []

    def __call__(self, *args, **kwargs):
        self.calls.append((args, kwargs))


class _StubInstance:
    """A class instance whose every accessed method is a fresh recorder.

    ``order`` is a shared list onto which each method name is appended in call
    order, so a wrapper's orchestration sequence can be asserted. Specific
    attributes (e.g. ``to_df``'s return value) can be pre-seeded by subclasses.
    """

    def __init__(self, order):
        object.__setattr__(self, "_order", order)
        object.__setattr__(self, "_recorders", {})

    def __getattr__(self, name):
        recorders = object.__getattribute__(self, "_recorders")
        order = object.__getattribute__(self, "_order")
        if name not in recorders:
            def _method(*args, _name=name, **kwargs):
                order.append(_name)
                return None
            recorders[name] = _method
        return recorders[name]


# ===========================================================================
# combine_feature_files() -- real merge logic, no mocking.
# ===========================================================================

# Columns dropped from df_ref before the merge on 'id' (from dicast.py:87).
_REF_DROP = ['sample', 'sv_type', 'chrom', 'chrom_2', 'start', 'end',
             'cohort', 'technology', 'caller', 'reference']
# Columns dropped from df_aln_ill before the merge on 'id' (from dicast.py:88).
_ALN_DROP = ['sample', 'sv_type', 'chrom', 'chrom_2', 'start', 'end', 'sv_len',
             'cohort', 'technology', 'caller', 'reference']


def _write_combine_inputs(tmp_path, sample, ref, ids):
    """Write the three feature TSVs combine_feature_files reads.

    Each file gets the drop-columns + 'id' + one unique kept column so we can
    assert exactly which columns survive the merge.
    """
    n = len(ids)
    workdir = tmp_path

    # raw TSV: 'id' + the kept-from-raw payload columns (nothing is dropped
    # from df_raw, so every column here survives).
    df_raw = pd.DataFrame({
        "id": ids,
        "raw_feature": list(range(n)),
        "sample": [sample] * n,
    })
    df_raw.to_csv(workdir / f"{sample}_{ref}.SVs.raw.tsv", sep="\t", index=False)

    # ref TSV: 'id' + every drop column + one unique kept column 'ref_feature'.
    ref_cols = {"id": ids, "ref_feature": [x * 10 for x in range(n)]}
    for c in _REF_DROP:
        ref_cols[c] = [f"r_{c}"] * n
    pd.DataFrame(ref_cols).to_csv(
        workdir / f"{sample}_{ref}.SVs.ref.tsv", sep="\t", index=False)

    # aln TSV (globbed): 'id' + every drop column + one unique kept column. Split
    # across two glob-matched files to exercise the pd.concat over the glob.
    half = max(1, n // 2)
    for part, sub_ids, sub_start, chrom in (
        (0, ids[:half], 0, "chr1"),
        (1, ids[half:], half, "chr2"),
    ):
        if not sub_ids:
            continue
        aln_cols = {
            "id": sub_ids,
            "aln_feature": [x * 100 for x in range(sub_start, sub_start + len(sub_ids))],
        }
        for c in _ALN_DROP:
            aln_cols[c] = [f"a_{c}"] * len(sub_ids)
        pd.DataFrame(aln_cols).to_csv(
            workdir / f"{sample}_{ref}.SVs.aln.ill.{chrom}.DEL.tsv",
            sep="\t", index=False)

    return str(workdir)


@pytest.mark.unit
def test_combine_feature_files_merges_on_id(tmp_path):
    sample, ref = "S1", "hg38"
    ids = ["v1", "v2", "v3", "v4"]
    workdir = _write_combine_inputs(tmp_path, sample, ref, ids)

    df = dicast.combine_feature_files(sample, ref, workdir)

    # Inner merge on 'id' with fully overlapping ids -> one row per id.
    assert len(df) == len(ids)
    assert sorted(df["id"]) == sorted(ids)

    # Kept payload columns from all three sources survive.
    for col in ("raw_feature", "ref_feature", "aln_feature"):
        assert col in df.columns

    # Drop columns coming from ref/aln must NOT appear (they were dropped before
    # the merge). 'sample' is special: it is dropped from BOTH ref and aln but
    # is a payload column of df_raw, so it DOES survive (from raw only).
    assert "sample" in df.columns
    for col in _REF_DROP:
        if col == "sample":
            continue
        assert col not in df.columns, f"{col!r} should have been dropped"

    # No duplicate / suffixed columns from the merge (e.g. 'sample_x').
    assert not any(c.endswith(("_x", "_y")) for c in df.columns)

    # Spot-check the join carried the right values together for one id.
    row = df[df["id"] == "v1"].iloc[0]
    assert row["raw_feature"] == 0
    assert row["ref_feature"] == 0
    assert row["aln_feature"] == 0


@pytest.mark.unit
def test_combine_feature_files_inner_join_drops_unmatched(tmp_path):
    # raw has an extra id that is absent from ref/aln -> inner join drops it.
    sample, ref = "S2", "hg38"
    ids = ["v1", "v2"]
    _write_combine_inputs(tmp_path, sample, ref, ids)

    # Append an extra raw-only id after the helper wrote the files.
    raw_path = tmp_path / f"{sample}_{ref}.SVs.raw.tsv"
    df_raw = pd.read_csv(raw_path, sep="\t")
    df_raw = pd.concat(
        [df_raw, pd.DataFrame({"id": ["orphan"], "raw_feature": [99],
                               "sample": [sample]})],
        ignore_index=True,
    )
    df_raw.to_csv(raw_path, sep="\t", index=False)

    df = dicast.combine_feature_files(sample, ref, str(tmp_path))
    assert "orphan" not in set(df["id"])
    assert sorted(df["id"]) == ["v1", "v2"]


@pytest.mark.unit
def test_combine_feature_files_preserves_string_sample_and_cohort_dtypes(tmp_path):
    """dicast.py:80 reads the raw TSV with dtype={'sample': str,
    'cohort_samples': str}. Without that, pandas' C parser would infer an
    all-digit column (a numeric-looking sample name, here '007') as int64 and
    silently drop the leading zero; a 'cohort_samples' column holding this
    repo's list-of-dicts repr string (see task brief) must also survive as a
    plain string rather than being touched by type inference."""
    sample, ref = "007", "hg38"
    ids = ["v1", "v2"]
    workdir = _write_combine_inputs(tmp_path, sample, ref, ids)

    raw_path = tmp_path / f"{sample}_{ref}.SVs.raw.tsv"
    # dtype=str here too: re-reading without it would itself re-typify '007'
    # as int64 before we even get to write the file combine_feature_files reads.
    df_raw = pd.read_csv(raw_path, sep="\t", dtype={"sample": str})
    cohort_repr = "[{'sample': '007', 'gt': '0/1'}]"
    df_raw["cohort_samples"] = cohort_repr
    df_raw.to_csv(raw_path, sep="\t", index=False)

    df = dicast.combine_feature_files(sample, ref, workdir)

    assert df["sample"].dtype == object
    assert (df["sample"] == "007").all()
    assert df["cohort_samples"].dtype == object
    assert (df["cohort_samples"] == cohort_repr).all()


# ===========================================================================
# add_info_tag_to_vcf() -- real vcfpy round-trip.
# ===========================================================================

_VCF_HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
"""


def _write_vcf(tmp_path, name, records):
    """Write a minimal vcfpy-parseable VCF. ``records`` is a list of
    (id, pos) tuples; each is a simple DEL record."""
    lines = [_VCF_HEADER.rstrip("\n")]
    for rec_id, pos in records:
        lines.append(
            f"chr1\t{pos}\t{rec_id}\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL\tGT\t0/1")
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\n")
    return str(path)


@pytest.mark.unit
def test_add_info_tag_to_vcf_writes_dq(tmp_path):
    sample, ref = "S1", "hg38"
    caller = "manta"

    # Two variants have scores; a third (UNSCORED) is intentionally absent from
    # the scores TSV to exercise the -1 fallback path.
    vcf_in = _write_vcf(
        tmp_path, "manta.vcf",
        [("DEL1", 1000), ("DEL2", 2000), ("UNSCORED", 3000)],
    )

    scores = pd.DataFrame({
        "id": ["DEL1", "DEL2"],
        "caller": [caller, caller],
        "dicast_qual": [0.91, 0.42],
        "sample": [sample, sample],
    })
    scores.to_csv(
        tmp_path / f"{sample}_{ref}.SVs.dicast.tsv", sep="\t", index=False)

    args = SimpleNamespace(
        workdir=str(tmp_path), sample=sample, ref=ref,
        vcfs=[(caller, vcf_in)],
    )

    dicast.add_info_tag_to_vcf(args)

    # Output lands in workdir, named after the input with '.vcf' -> '.dicast.vcf'
    # (workdir == tmp_path here, where the input also lives).
    out_path = vcf_in.replace(".vcf", ".dicast.vcf")
    reader = vcfpy.Reader.from_path(out_path)

    # The DQ INFO line must have been added to the header.
    assert "DQ" in reader.header.info_ids()

    dq_by_id = {}
    for rec in reader:
        dq_by_id[rec.ID[0]] = rec.INFO["DQ"]

    # Scored variants carry their dicast_qual as a string; the unscored one
    # falls back to '-1'.
    assert dq_by_id["DEL1"] == "0.91"
    assert dq_by_id["DEL2"] == "0.42"
    assert dq_by_id["UNSCORED"] == "-1"


@pytest.mark.unit
def test_add_info_tag_to_vcf_filters_by_caller(tmp_path):
    # A scores TSV row for a different caller must NOT match this VCF's records,
    # so every record falls back to '-1'.
    sample, ref = "S1", "hg38"
    vcf_in = _write_vcf(tmp_path, "delly.vcf", [("DEL1", 1000)])

    pd.DataFrame({
        "id": ["DEL1"],
        "caller": ["manta"],          # different caller than 'delly'
        "dicast_qual": [0.99],
        "sample": [sample],
    }).to_csv(tmp_path / f"{sample}_{ref}.SVs.dicast.tsv", sep="\t", index=False)

    args = SimpleNamespace(
        workdir=str(tmp_path), sample=sample, ref=ref,
        vcfs=[("delly", vcf_in)],
    )
    dicast.add_info_tag_to_vcf(args)

    reader = vcfpy.Reader.from_path(vcf_in.replace(".vcf", ".dicast.vcf"))
    rec = next(iter(reader))
    assert rec.INFO["DQ"] == "-1"


# ===========================================================================
# extract_reference_features() -- heavy ReferenceAnnotator mocked.
# ===========================================================================

@pytest.mark.unit
def test_extract_reference_features_orchestration(tmp_path, monkeypatch):
    order = []
    constructed = {}

    class _FakeRA(_StubInstance):
        def __init__(self, reference_filenames):
            super().__init__(order)
            constructed["reference_filenames"] = reference_filenames
            # extract_reference_features prints RA.df_calls_annot.shape twice.
            object.__setattr__(
                self, "df_calls_annot",
                SimpleNamespace(shape=(3, 4)))

    monkeypatch.setattr(dicast, "ReferenceAnnotator", _FakeRA)

    args = SimpleNamespace(
        workdir=str(tmp_path), ref="hg38",
        repeats="/r/repeats", vntrs="/r/vntrs", strs="/r/strs",
        cgis="/r/cgis", centromeres="/r/centromeres", gaps="/r/gaps",
        althaps="/r/althaps", gc="/r/gc",
    )

    dicast.extract_reference_features(args, "S1")

    # The reference filename dict is assembled from the args (verify the mapping
    # is wired correctly for a couple of keys).
    assert constructed["reference_filenames"]["repeats_filename"] == "/r/repeats"
    assert constructed["reference_filenames"]["gc_filename"] == "/r/gc"
    assert constructed["reference_filenames"]["cpgislands_filename"] == "/r/cgis"

    # The full annotation pipeline runs in this exact order (dicast.py:147-169).
    assert order == [
        "load_from_csv",
        "split_bnd",
        "annotate_repeats",
        "annotate_vntrs",
        "annotate_strs",
        "annotate_cpg_islands",
        "annotate_centromeres",
        "annotate_asmb_gaps",
        "annotate_alt_haps",
        "annotate_gc_content",
        "aggregate_results",
        "to_csv",
    ]


# ===========================================================================
# collect_aln_features() -- heavy AlignmentAnnotatorIllumina mocked.
# ===========================================================================

@pytest.mark.unit
def test_collect_aln_features_orchestration(monkeypatch):
    order = []
    constructed = {}

    class _FakeAAI(_StubInstance):
        def __init__(self, bam_filename, chrom, sv_type, sample):
            super().__init__(order)
            constructed.update(
                bam_filename=bam_filename, chrom=chrom,
                sv_type=sv_type, sample=sample)

    monkeypatch.setattr(dicast, "AlignmentAnnotatorIllumina", _FakeAAI)

    dicast.collect_aln_features(
        bam_filename="/x/sample.bam",
        variant_filename="/w/in.tsv",
        variant_annot_filename="/w/out.tsv",
        chrom="chr1",
        sv_type="DEL",
        sample="S1",
    )

    # Constructor positional args are forwarded straight through.
    assert constructed == {
        "bam_filename": "/x/sample.bam", "chrom": "chr1",
        "sv_type": "DEL", "sample": "S1",
    }

    # Methods invoked in the documented order (dicast.py:184-191).
    assert order == [
        "load_from_csv",
        "calculate_coverage_baseline",
        "calculate_insertsize_baseline",
        "calculate_mapping_quality_baseline",
        "annotate_coverage",
        "annotate_read_based_features",
        "to_csv",
    ]


# ===========================================================================
# score_variants() -- heavy Dicast model mocked.
# ===========================================================================

def _make_fake_dicast(per_type_orders, loaded_models):
    """Builds a fake Dicast class recording, per sv_type, its ordered method
    calls (into ``per_type_orders``) and every model filename passed to
    ``load()`` (into ``loaded_models``).

    Note: no ``score_translocations`` method is defined here. In this repo
    ``score_variants`` only special-cases 'INV' -- every other sv_type,
    'BND' included, goes through the model-loading branch below (see module
    docstring's "Dropped" section).
    """
    class _FakeDicast:
        def __init__(self, sv_type):
            self.sv_type = sv_type
            self._order = per_type_orders.setdefault(sv_type, [])

        def load_from_csv(self, fn):
            self._order.append("load_from_csv")

        def impute_missing_values(self):
            self._order.append("impute_missing_values")

        def load(self, model_filename):
            self._order.append("load")
            loaded_models.append(model_filename)

        def predict(self):
            self._order.append("predict")

        def score_inversions(self):
            self._order.append("score_inversions")

        def to_df(self):
            self._order.append("to_df")
            # One row tagged with the sv_type so we can verify concat output.
            return pd.DataFrame({
                "id": [f"{self.sv_type}_1"],
                "sv_type": [self.sv_type],
                "dicast_qual": [0.5],
            })
    return _FakeDicast


@pytest.mark.unit
def test_score_variants_dispatches_per_sv_type(tmp_path, monkeypatch):
    # Record, per sv_type, the ordered method calls and the model filename used.
    per_type_orders = {}
    loaded_models = []
    monkeypatch.setattr(dicast, "Dicast", _make_fake_dicast(per_type_orders, loaded_models))

    # args has no 'pop' attribute at all -- getattr(arguments, 'pop', False)
    # must fall back to False (dicast.py:204) for backward compatibility.
    args = SimpleNamespace(
        workdir=str(tmp_path), ref="hg38", models="/models",
    )
    sv_types = ["DEL", "DUP", "INV"]

    dicast.score_variants(sv_types, args, "S1")

    # DEL and DUP are XGBoost model types: each loads its own model JSON and
    # predicts.
    assert per_type_orders["DEL"] == [
        "load_from_csv", "impute_missing_values", "load", "predict", "to_df"]
    assert per_type_orders["DUP"] == [
        "load_from_csv", "impute_missing_values", "load", "predict", "to_df"]
    assert loaded_models == ["/models/dicast_DEL.json", "/models/dicast_DUP.json"]

    # INV is scored heuristically (no model load).
    assert per_type_orders["INV"] == [
        "load_from_csv", "impute_missing_values", "score_inversions", "to_df"]

    # The concatenated predictions are written to the dicast TSV.
    out = pd.read_csv(tmp_path / "S1_hg38.SVs.dicast.tsv", sep="\t")
    assert sorted(out["sv_type"]) == ["DEL", "DUP", "INV"]
    assert len(out) == 3


@pytest.mark.unit
def test_score_variants_pop_flag_uses_pop_model_when_present(tmp_path, monkeypatch):
    per_type_orders = {}
    loaded_models = []
    monkeypatch.setattr(dicast, "Dicast", _make_fake_dicast(per_type_orders, loaded_models))

    models_dir = tmp_path / "models"
    models_dir.mkdir()
    (models_dir / "dicast_DEL_pop.json").write_text("{}")

    args = SimpleNamespace(
        workdir=str(tmp_path), ref="hg38", models=str(models_dir), pop=True,
    )
    dicast.score_variants(["DEL"], args, "S1")

    assert loaded_models == [str(models_dir / "dicast_DEL_pop.json")]


@pytest.mark.unit
def test_score_variants_pop_flag_falls_back_when_pop_model_missing(tmp_path, monkeypatch):
    per_type_orders = {}
    loaded_models = []
    monkeypatch.setattr(dicast, "Dicast", _make_fake_dicast(per_type_orders, loaded_models))

    models_dir = tmp_path / "models"
    models_dir.mkdir()
    # No dicast_INS_pop.json written -> must fall back to the standard model.

    args = SimpleNamespace(
        workdir=str(tmp_path), ref="hg38", models=str(models_dir), pop=True,
    )
    dicast.score_variants(["INS"], args, "S1")

    assert loaded_models == [str(models_dir / "dicast_INS.json")]


@pytest.mark.unit
def test_score_variants_pop_flag_only_applies_to_del_and_ins(tmp_path, monkeypatch):
    # Even with a pop model file present for DUP, --pop must not pick it up:
    # dicast.py:209 only special-cases sv_type in ('DEL', 'INS').
    per_type_orders = {}
    loaded_models = []
    monkeypatch.setattr(dicast, "Dicast", _make_fake_dicast(per_type_orders, loaded_models))

    models_dir = tmp_path / "models"
    models_dir.mkdir()
    (models_dir / "dicast_DUP_pop.json").write_text("{}")

    args = SimpleNamespace(
        workdir=str(tmp_path), ref="hg38", models=str(models_dir), pop=True,
    )
    dicast.score_variants(["DUP"], args, "S1")

    assert loaded_models == [str(models_dir / "dicast_DUP.json")]


# ===========================================================================
# Integration: the real `call` pipeline end-to-end on the shipped demo data.
# ===========================================================================

TEST_DATA_DIR = REPO_DIR / "tests/data"


@pytest.mark.integration
def test_call_pipeline_end_to_end_on_demo_data(tmp_path):
    """Runs `python3 dicast.py call` on the shipped chr21 demo dataset, with
    nothing mocked. Exercises the whole pipeline (VariantPrep -> reference
    feature extraction -> alignment feature extraction -> combine ->
    score_variants -> add_info_tag_to_vcf) against real code.

    Runtime is a few seconds -- the demo BAM is intentionally downsampled at
    each of the 20 variant loci to ~1-2x specifically so no call is dropped by
    the coverage-outlier filter (cov_thr=5); see tests/data/README.md's "Why
    the variant loci are downsampled this hard" for the full explanation.

    All outputs, including the DQ-tagged `*.dicast.vcf`, land in `--workdir`,
    so the repo's own tests/data/ inputs can be used in place.
    """
    workdir = tmp_path / "workdir"
    workdir.mkdir()

    vcf_in = TEST_DATA_DIR / "demo_delly.vcf.gz"

    cmd = [
        sys.executable, "-m", "dicast", "call",
        "--sample", "demo",
        "--workdir", str(workdir),
        "--fai", str(TEST_DATA_DIR / "hg38.fa.fai"),
        "--bam", str(TEST_DATA_DIR / "demo.bam"),
        "--vcfs", f"delly={vcf_in}",
        "--annot-dir", str(TEST_DATA_DIR / "annot"),
        "--chrom", "chr21",
        "--sv_types", "DEL", "INS",
        "--threads", "2",
    ]
    result = subprocess.run(
        cmd, cwd=str(REPO_DIR), capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, (
        f"dicast.py call exited {result.returncode}\n"
        f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}")

    scores_path = workdir / "demo_hg38.SVs.dicast.tsv"
    assert scores_path.is_file()
    df = pd.read_csv(scores_path, sep="\t")

    # All 20 shipped calls (10 DEL, 10 INS) survive feature collection and
    # scoring -- none get dropped as a coverage-outlier artifact.
    assert len(df) == 20
    assert df["id"].is_unique
    assert sorted(df["sv_type"].unique()) == ["DEL", "INS"]
    assert set(df["caller"]) == {"delly"}

    # Every call got a real (non-fallback) score in [0, 1] from the shipped
    # models, not a missing/NA value.
    assert df["dicast_qual"].notna().all()
    assert df["dicast_qual"].between(0, 1).all()

    # The VCF got its DQ tags too, matching the scores TSV; the re-emitted
    # VCF is written into --workdir, not next to the input.
    dq_vcf_path = str(workdir / "demo_delly.dicast.vcf")
    reader = vcfpy.Reader.from_path(dq_vcf_path)
    scores_by_id = dict(zip(df["id"], df["dicast_qual"]))
    seen = 0
    for rec in reader:
        assert rec.ID[0] in scores_by_id
        assert float(rec.INFO["DQ"]) == pytest.approx(scores_by_id[rec.ID[0]])
        seen += 1
    assert seen == 20
