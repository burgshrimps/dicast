"""Unit tests for ``dicast/parsing.py``.

This repo's CLI differs substantially from the older lucid/dicast dev line this
module was ported from:

- There is no ``cohort`` subcommand any more. Cohort/CSV cohort mode was
  replaced by a ``multi`` subcommand (cross-sample rescue, e.g. a trio) whose
  ``--bams``/``--vcfs`` tokens have their own ``sample=bam_file`` /
  ``sample:caller=vcf_file`` formats.
- There is no exome support at all (no ``--exome``/``--exome_regions`` flags).
- ``call`` gained many required/new flags: ``--sample``, ``--bam`` and
  ``--vcfs`` are now ``required=True``; ``--annot-dir`` resolves individual
  annotation-file flags to canonical hg38 filenames unless explicitly
  overridden; ``--models`` has a package-relative default instead of being
  required; and ``--pop``/``--pop-catalog``/``--benchmark``/``--sv_types`` are
  new.
- ``parse_arguments`` now also runs ``resolve_annotation_paths`` on the parsed
  namespace before returning it, and a separate ``validate_inputs`` function
  (not exercised by ``parse_arguments`` itself) does on-disk validation
  (missing files, missing BAM indexes, missing model files, malformed
  tokens, ...).
- There is a third subcommand, ``check`` (report what dicast understands
  from a set of input VCFs, without running the pipeline): ``--vcfs`` and
  ``--fai`` are required, ``--sample``/``--chrom``/``--out`` are optional.
  It has no workdir/models/annotation flags at all, so
  ``resolve_annotation_paths`` is a no-op for it (it returns early for any
  command other than 'call'/'multi') and its own ``_validate_check_inputs``
  only checks that ``--fai`` and every ``--vcfs`` file exist.
- ``parse_arguments`` is now ``build_parser()`` (returns the configured
  ``argparse.ArgumentParser``, all three subcommands) plus a thin
  ``parse_arguments()`` that parses argv against it and resolves annotation
  paths.

Expected values are derived from the argv list we construct, or from real
files created in ``tmp_path``, not from the function's own output.
"""
from __future__ import annotations

import argparse
import os

import pytest

from dicast import parsing


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _touch(path):
    """Create an empty file (and parent dirs) at ``path`` and return its str."""
    path = str(path)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w'):
        pass
    return path


@pytest.fixture
def call_required_args(tmp_path, monkeypatch):
    """Minimal argv tokens satisfying every required flag of the 'call' subcommand.

    Uses an empty, isolated --annot-dir/--models so tests don't depend on (or
    slow down against) the real multi-gigabyte annotation files/model files.
    DICAST_DATA_DIR points at a pre-populated tmp store so tests that omit
    --annot-dir resolve against it and never trigger the first-use download.
    """
    fai = _touch(tmp_path / "ref.fa.fai")
    bam = _touch(tmp_path / "sample.bam")
    _touch(tmp_path / "sample.bam.bai")
    vcf = _touch(tmp_path / "manta.vcf")
    annot_dir = tmp_path / "annot_empty"
    annot_dir.mkdir()
    models_dir = tmp_path / "models_empty"
    models_dir.mkdir()
    data_dir = tmp_path / "data_store"
    data_dir.mkdir()
    _populate_annot_dir(data_dir)
    monkeypatch.setenv("DICAST_DATA_DIR", str(data_dir))
    return {
        "fai": fai,
        "bam": bam,
        "vcf": vcf,
        "annot_dir": str(annot_dir),
        "models_dir": str(models_dir),
        "data_dir": str(data_dir),
        "workdir": str(tmp_path / "work"),
    }


def _populate_annot_dir(annot_dir):
    """Create empty files under annot_dir for every canonical annotation name."""
    for filename in parsing.ANNOT_CANONICAL_NAMES.values():
        _touch(os.path.join(annot_dir, filename))


def _populate_models_dir(models_dir, sv_types=('DEL', 'DUP', 'INS')):
    for sv_type in sv_types:
        _touch(os.path.join(models_dir, f'dicast_{sv_type}.json'))


# ---------------------------------------------------------------------------
# call subcommand — valid invocations
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_call_subcommand_parses_fields(call_required_args):
    args = parsing.parse_arguments(
        arguments=[
            "call",
            "--cohort", "mycohort",
            "--sample", "mysample",
            "--chrom", "chr1", "chr2",
            "--ref", "hg38",
            "--technology", "ill",
            "--workdir", call_required_args["workdir"],
            "--fai", call_required_args["fai"],
            "--bam", call_required_args["bam"],
            "--vcfs", f"manta={call_required_args['vcf']}",
            "--threads", "4",
        ]
    )
    assert args.command == "call"
    assert args.cohort == "mycohort"
    assert args.sample == "mysample"
    # --chrom uses nargs='+', so it collects into a list.
    assert args.chrom == ["chr1", "chr2"]
    assert args.ref == "hg38"
    assert args.technology == "ill"
    assert args.workdir == call_required_args["workdir"]
    assert args.fai == call_required_args["fai"]
    assert args.bam == call_required_args["bam"]
    # --vcfs uses type=lambda kv: kv.split('='), so each token becomes [method, file].
    assert args.vcfs == [["manta", call_required_args["vcf"]]]
    # --threads has type=int, so the value is coerced from the string "4".
    assert args.threads == 4
    assert isinstance(args.threads, int)


@pytest.mark.unit
def test_call_subcommand_defaults(call_required_args):
    # Supply only the required flags; every optional arg should fall back to
    # its declared default.
    args = parsing.parse_arguments(
        arguments=[
            "call",
            "--sample", "mysample",
            "--workdir", call_required_args["workdir"],
            "--fai", call_required_args["fai"],
            "--bam", call_required_args["bam"],
            "--vcfs", f"manta={call_required_args['vcf']}",
        ]
    )
    assert args.command == "call"
    # Declared defaults from the source.
    assert args.cohort == "none"
    assert args.chrom == "all"
    assert args.ref == "hg38"
    assert args.technology == "ill"
    assert args.threads == 1
    assert args.pop is False
    # --annot-dir was not passed, so it falls back to the managed annotation
    # store (isolated to a populated tmp dir via DICAST_DATA_DIR by the
    # fixture), and --pop-catalog is resolved against it like the other flags.
    assert args.annot_dir == call_required_args["data_dir"]
    assert args.pop_catalog == os.path.join(args.annot_dir, parsing.POP_CATALOG_NAME)
    assert args.benchmark is None
    assert args.sv_types is None
    # --models defaults to the package's own models directory.
    assert args.models == parsing.DEFAULT_MODELS_DIR
    # This CLI has no exome support at all.
    assert not hasattr(args, 'exome')
    assert not hasattr(args, 'exome_regions')


@pytest.mark.unit
def test_call_vcfs_split_on_equals(call_required_args):
    # --vcfs uses type=lambda kv: kv.split('=') so each token becomes [method, file].
    vcf2 = _touch(os.path.dirname(call_required_args["vcf"]) + "/delly.vcf")
    args = parsing.parse_arguments(
        arguments=[
            "call",
            "--sample", "mysample",
            "--workdir", call_required_args["workdir"],
            "--fai", call_required_args["fai"],
            "--bam", call_required_args["bam"],
            "--vcfs", f"manta={call_required_args['vcf']}", f"delly={vcf2}",
        ]
    )
    assert args.vcfs == [["manta", call_required_args["vcf"]], ["delly", vcf2]]


# ---------------------------------------------------------------------------
# call subcommand — required arguments (new in this CLI)
# ---------------------------------------------------------------------------

@pytest.mark.unit
@pytest.mark.parametrize("missing", ["--sample", "--workdir", "--fai", "--bam", "--vcfs"])
def test_call_missing_required_flag_exits(call_required_args, missing):
    # --sample/--workdir/--fai/--bam/--vcfs are all required=True in this CLI
    # (unlike the older dev line, where only --sample/--workdir/--fai were
    # required and there was no --bam/--vcfs requirement check at parse time).
    full_args = [
        "call",
        "--sample", "mysample",
        "--workdir", call_required_args["workdir"],
        "--fai", call_required_args["fai"],
        "--bam", call_required_args["bam"],
        "--vcfs", f"manta={call_required_args['vcf']}",
    ]
    # Drop the flag (and its value) under test.
    idx = full_args.index(missing)
    del full_args[idx:idx + 2]
    with pytest.raises(SystemExit):
        parsing.parse_arguments(arguments=full_args)


# ---------------------------------------------------------------------------
# --annot-dir canonical-filename resolution / explicit-flag override
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_annot_dir_fills_in_canonical_paths(call_required_args, tmp_path):
    annot_dir = call_required_args["annot_dir"]
    args = parsing.parse_arguments(
        arguments=[
            "call",
            "--sample", "mysample",
            "--workdir", call_required_args["workdir"],
            "--fai", call_required_args["fai"],
            "--bam", call_required_args["bam"],
            "--vcfs", f"manta={call_required_args['vcf']}",
            "--annot-dir", annot_dir,
        ]
    )
    for flag, filename in parsing.ANNOT_CANONICAL_NAMES.items():
        assert getattr(args, flag) == os.path.join(annot_dir, filename)
    assert args.pop_catalog == os.path.join(annot_dir, parsing.POP_CATALOG_NAME)


@pytest.mark.unit
def test_explicit_annot_flag_overrides_annot_dir(call_required_args, tmp_path):
    annot_dir = call_required_args["annot_dir"]
    custom_repeats = _touch(tmp_path / "custom_repeats.tsv")
    args = parsing.parse_arguments(
        arguments=[
            "call",
            "--sample", "mysample",
            "--workdir", call_required_args["workdir"],
            "--fai", call_required_args["fai"],
            "--bam", call_required_args["bam"],
            "--vcfs", f"manta={call_required_args['vcf']}",
            "--annot-dir", annot_dir,
            "--repeats", custom_repeats,
        ]
    )
    # The explicitly passed flag wins over the --annot-dir-derived default.
    assert args.repeats == custom_repeats
    # Untouched flags are still filled in from --annot-dir.
    assert args.cgis == os.path.join(annot_dir, parsing.ANNOT_CANONICAL_NAMES['cgis'])


# ---------------------------------------------------------------------------
# --pop / --pop-catalog
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_pop_catalog_defaults_from_annot_dir(call_required_args):
    annot_dir = call_required_args["annot_dir"]
    args = parsing.parse_arguments(
        arguments=[
            "call",
            "--sample", "mysample",
            "--workdir", call_required_args["workdir"],
            "--fai", call_required_args["fai"],
            "--bam", call_required_args["bam"],
            "--vcfs", f"manta={call_required_args['vcf']}",
            "--annot-dir", annot_dir,
        ]
    )
    assert args.pop_catalog == os.path.join(annot_dir, parsing.POP_CATALOG_NAME)


@pytest.mark.unit
def test_pop_flag_appends_pav_catalog_to_call_vcfs(call_required_args):
    annot_dir = call_required_args["annot_dir"]
    args = parsing.parse_arguments(
        arguments=[
            "call",
            "--sample", "mysample",
            "--workdir", call_required_args["workdir"],
            "--fai", call_required_args["fai"],
            "--bam", call_required_args["bam"],
            "--vcfs", f"manta={call_required_args['vcf']}",
            "--annot-dir", annot_dir,
            "--pop",
        ]
    )
    assert args.pop is True
    expected_catalog = os.path.join(annot_dir, parsing.POP_CATALOG_NAME)
    assert args.vcfs == [
        ["manta", call_required_args["vcf"]],
        ["pav", expected_catalog],
    ]


@pytest.mark.unit
def test_pop_flag_uses_explicit_pop_catalog(call_required_args):
    custom_catalog = os.path.dirname(call_required_args["vcf"]) + "/custom_pav.vcf.gz"
    _touch(custom_catalog)
    args = parsing.parse_arguments(
        arguments=[
            "call",
            "--sample", "mysample",
            "--workdir", call_required_args["workdir"],
            "--fai", call_required_args["fai"],
            "--bam", call_required_args["bam"],
            "--vcfs", f"manta={call_required_args['vcf']}",
            "--pop",
            "--pop-catalog", custom_catalog,
        ]
    )
    assert args.pop_catalog == custom_catalog
    assert ["pav", custom_catalog] in args.vcfs


# ---------------------------------------------------------------------------
# --models default
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_models_default_is_repo_relative(call_required_args):
    args = parsing.parse_arguments(
        arguments=[
            "call",
            "--sample", "mysample",
            "--workdir", call_required_args["workdir"],
            "--fai", call_required_args["fai"],
            "--bam", call_required_args["bam"],
            "--vcfs", f"manta={call_required_args['vcf']}",
        ]
    )
    assert args.models == parsing.DEFAULT_MODELS_DIR


@pytest.mark.unit
def test_models_explicit_flag_overrides_default(call_required_args):
    args = parsing.parse_arguments(
        arguments=[
            "call",
            "--sample", "mysample",
            "--workdir", call_required_args["workdir"],
            "--fai", call_required_args["fai"],
            "--bam", call_required_args["bam"],
            "--vcfs", f"manta={call_required_args['vcf']}",
            "--models", call_required_args["models_dir"],
        ]
    )
    assert args.models == call_required_args["models_dir"]


# ---------------------------------------------------------------------------
# validate_inputs — 'call' subcommand
# ---------------------------------------------------------------------------

def _build_valid_call_args(call_required_args, tmp_path):
    _populate_annot_dir(call_required_args["annot_dir"])
    _populate_models_dir(call_required_args["models_dir"], sv_types=('DEL', 'DUP', 'INS'))
    return parsing.parse_arguments(
        arguments=[
            "call",
            "--sample", "mysample",
            "--workdir", call_required_args["workdir"],
            "--fai", call_required_args["fai"],
            "--bam", call_required_args["bam"],
            "--vcfs", f"manta={call_required_args['vcf']}",
            "--annot-dir", call_required_args["annot_dir"],
            "--models", call_required_args["models_dir"],
        ]
    )


@pytest.mark.unit
def test_validate_call_inputs_all_present_succeeds(call_required_args, tmp_path):
    args = _build_valid_call_args(call_required_args, tmp_path)
    # Should not raise/exit, and should create the workdir.
    parsing.validate_inputs(args)
    assert os.path.isdir(args.workdir)


@pytest.mark.unit
def test_validate_call_inputs_missing_bam_exits(call_required_args, tmp_path, capsys):
    args = _build_valid_call_args(call_required_args, tmp_path)
    os.remove(args.bam)
    with pytest.raises(SystemExit):
        parsing.validate_inputs(args)
    err = capsys.readouterr().err
    assert '--bam file' in err
    assert args.bam in err


@pytest.mark.unit
def test_validate_call_inputs_missing_fai_exits(call_required_args, tmp_path, capsys):
    args = _build_valid_call_args(call_required_args, tmp_path)
    os.remove(args.fai)
    with pytest.raises(SystemExit):
        parsing.validate_inputs(args)
    err = capsys.readouterr().err
    assert '--fai file' in err


@pytest.mark.unit
def test_validate_call_inputs_missing_annot_file_exits(call_required_args, tmp_path, capsys):
    args = _build_valid_call_args(call_required_args, tmp_path)
    os.remove(args.cgis)
    with pytest.raises(SystemExit):
        parsing.validate_inputs(args)
    err = capsys.readouterr().err
    assert '--cgis file not found' in err


@pytest.mark.unit
def test_validate_call_inputs_missing_vcf_exits(call_required_args, tmp_path, capsys):
    args = _build_valid_call_args(call_required_args, tmp_path)
    os.remove(call_required_args["vcf"])
    with pytest.raises(SystemExit):
        parsing.validate_inputs(args)
    err = capsys.readouterr().err
    assert 'VCF file for caller manta not found' in err


@pytest.mark.unit
def test_validate_call_inputs_no_vcfs_reports_required(call_required_args, tmp_path, capsys):
    args = _build_valid_call_args(call_required_args, tmp_path)
    args.vcfs = []
    with pytest.raises(SystemExit):
        parsing.validate_inputs(args)
    err = capsys.readouterr().err
    assert '--vcfs is required' in err


@pytest.mark.unit
def test_validate_call_inputs_malformed_vcfs_token_exits(call_required_args, tmp_path, capsys):
    # A --vcfs token without '=' (e.g. "mantavcf") splits into a single-element
    # list, which the validator flags as malformed rather than crashing.
    args = _build_valid_call_args(call_required_args, tmp_path)
    args.vcfs = [["mantavcf"]]
    with pytest.raises(SystemExit):
        parsing.validate_inputs(args)
    err = capsys.readouterr().err
    assert 'Malformed --vcfs entry' in err
    assert 'mantavcf' in err


@pytest.mark.unit
def test_validate_call_inputs_missing_bam_index_exits(call_required_args, tmp_path, capsys):
    args = _build_valid_call_args(call_required_args, tmp_path)
    os.remove(args.bam + '.bai')
    with pytest.raises(SystemExit):
        parsing.validate_inputs(args)
    err = capsys.readouterr().err
    assert 'BAM index not found' in err
    assert 'samtools index' in err


@pytest.mark.unit
def test_validate_call_inputs_missing_model_file_exits(call_required_args, tmp_path, capsys):
    args = _build_valid_call_args(call_required_args, tmp_path)
    os.remove(os.path.join(args.models, 'dicast_DUP.json'))
    with pytest.raises(SystemExit):
        parsing.validate_inputs(args)
    err = capsys.readouterr().err
    assert 'Model file not found for SV type DUP' in err


@pytest.mark.unit
def test_validate_call_inputs_missing_annot_files_hint_fetch_command(call_required_args, tmp_path, capsys):
    # Missing annotation files hint at the dicast-fetch-annotations command.
    args = _build_valid_call_args(call_required_args, tmp_path)
    os.remove(args.gc)
    os.remove(args.repeats)
    with pytest.raises(SystemExit):
        parsing.validate_inputs(args)
    err = capsys.readouterr().err
    assert '--gc file not found' in err
    assert '--repeats file not found' in err
    assert err.count('dicast-fetch-annotations') == 2
    assert '--cgis file not found' not in err


@pytest.mark.unit
def test_validate_call_inputs_pop_catalog_missing_hints_fetch(call_required_args, tmp_path, capsys):
    args = _build_valid_call_args(call_required_args, tmp_path)
    args.pop = True
    args.pop_catalog = str(tmp_path / "missing_pav_catalog.vcf.gz")
    with pytest.raises(SystemExit):
        parsing.validate_inputs(args)
    err = capsys.readouterr().err
    assert '--pop-catalog file not found' in err
    assert 'dicast-fetch-annotations' in err


# ---------------------------------------------------------------------------
# multi subcommand — token parsing
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_multi_bams_split_on_equals(call_required_args, tmp_path):
    bam2 = _touch(tmp_path / "sample2.bam")
    args = parsing.parse_arguments(
        arguments=[
            "multi",
            "--workdir", call_required_args["workdir"],
            "--fai", call_required_args["fai"],
            "--bams", f"kid={call_required_args['bam']}", f"mom={bam2}",
            "--vcfs", f"kid:manta={call_required_args['vcf']}", f"mom:manta={call_required_args['vcf']}",
        ]
    )
    assert args.command == "multi"
    assert args.bams == [["kid", call_required_args["bam"]], ["mom", bam2]]


@pytest.mark.unit
def test_multi_vcfs_split_on_colon_and_equals(call_required_args, tmp_path):
    # --vcfs in 'multi' mode uses _parse_multi_vcf_token: 'sample:caller=vcf_file'.
    args = parsing.parse_arguments(
        arguments=[
            "multi",
            "--workdir", call_required_args["workdir"],
            "--fai", call_required_args["fai"],
            "--bams", f"kid={call_required_args['bam']}",
            "--vcfs", f"kid:manta={call_required_args['vcf']}",
        ]
    )
    assert args.vcfs == [["kid", "manta", call_required_args["vcf"]]]


@pytest.mark.unit
def test_multi_missing_required_flags_exits(call_required_args):
    with pytest.raises(SystemExit):
        parsing.parse_arguments(
            arguments=[
                "multi",
                "--workdir", call_required_args["workdir"],
                "--fai", call_required_args["fai"],
            ]
        )


@pytest.mark.unit
def test_multi_pop_expands_per_bam_sample(call_required_args, tmp_path):
    annot_dir = call_required_args["annot_dir"]
    bam2 = _touch(tmp_path / "sample2.bam")
    args = parsing.parse_arguments(
        arguments=[
            "multi",
            "--workdir", call_required_args["workdir"],
            "--fai", call_required_args["fai"],
            "--bams", f"kid={call_required_args['bam']}", f"mom={bam2}",
            "--vcfs", f"kid:manta={call_required_args['vcf']}", f"mom:manta={call_required_args['vcf']}",
            "--annot-dir", annot_dir,
            "--pop",
        ]
    )
    expected_catalog = os.path.join(annot_dir, parsing.POP_CATALOG_NAME)
    pav_entries = [entry for entry in args.vcfs if entry[1] == 'pav']
    assert sorted(e[0] for e in pav_entries) == ['kid', 'mom']
    assert all(e[2] == expected_catalog for e in pav_entries)


# ---------------------------------------------------------------------------
# validate_inputs — 'multi' subcommand
# ---------------------------------------------------------------------------

def _build_valid_multi_args(call_required_args, tmp_path, second_sample="mom"):
    _populate_annot_dir(call_required_args["annot_dir"])
    _populate_models_dir(call_required_args["models_dir"], sv_types=('DEL', 'DUP', 'INS'))
    bam2 = _touch(tmp_path / "sample2.bam")
    _touch(str(bam2) + '.bai')
    return parsing.parse_arguments(
        arguments=[
            "multi",
            "--workdir", call_required_args["workdir"],
            "--fai", call_required_args["fai"],
            "--bams", f"kid={call_required_args['bam']}", f"{second_sample}={bam2}",
            "--vcfs", f"kid:manta={call_required_args['vcf']}", f"{second_sample}:manta={call_required_args['vcf']}",
            "--annot-dir", call_required_args["annot_dir"],
            "--models", call_required_args["models_dir"],
        ]
    )


@pytest.mark.unit
def test_validate_multi_inputs_all_present_succeeds(call_required_args, tmp_path):
    args = _build_valid_multi_args(call_required_args, tmp_path)
    parsing.validate_inputs(args)
    assert os.path.isdir(args.workdir)


@pytest.mark.unit
def test_validate_multi_inputs_too_few_samples_exits(call_required_args, tmp_path, capsys):
    _populate_annot_dir(call_required_args["annot_dir"])
    _populate_models_dir(call_required_args["models_dir"], sv_types=('DEL', 'DUP', 'INS'))
    args = parsing.parse_arguments(
        arguments=[
            "multi",
            "--workdir", call_required_args["workdir"],
            "--fai", call_required_args["fai"],
            "--bams", f"kid={call_required_args['bam']}",
            "--vcfs", f"kid:manta={call_required_args['vcf']}",
            "--annot-dir", call_required_args["annot_dir"],
            "--models", call_required_args["models_dir"],
        ]
    )
    with pytest.raises(SystemExit):
        parsing.validate_inputs(args)
    err = capsys.readouterr().err
    assert 'at least two samples' in err


@pytest.mark.unit
def test_validate_multi_inputs_bams_vcfs_sample_mismatch_exits(call_required_args, tmp_path, capsys):
    # A --bams sample with no matching --vcfs entry, and a --vcfs sample with
    # no matching --bams entry, must both be reported.
    args = _build_valid_multi_args(call_required_args, tmp_path, second_sample="mom")
    # Rename the second --vcfs sample so it no longer matches any --bams sample.
    args.vcfs = [
        entry if entry[0] != 'mom' else ['dad', entry[1], entry[2]]
        for entry in args.vcfs
    ]
    with pytest.raises(SystemExit):
        parsing.validate_inputs(args)
    err = capsys.readouterr().err
    assert 'No --vcfs entries for sample(s) with a --bams entry: mom' in err
    assert '--vcfs references sample(s) without a --bams entry: dad' in err


@pytest.mark.unit
def test_validate_multi_inputs_malformed_bams_token_exits(call_required_args, tmp_path, capsys):
    args = _build_valid_multi_args(call_required_args, tmp_path)
    args.bams = [["onlysample"]]
    with pytest.raises(SystemExit):
        parsing.validate_inputs(args)
    err = capsys.readouterr().err
    assert 'Malformed --bams entry' in err


@pytest.mark.unit
def test_validate_multi_inputs_missing_bam_index_exits(call_required_args, tmp_path, capsys):
    args = _build_valid_multi_args(call_required_args, tmp_path)
    kid_bam = call_required_args["bam"]
    os.remove(kid_bam + '.bai')
    with pytest.raises(SystemExit):
        parsing.validate_inputs(args)
    err = capsys.readouterr().err
    assert 'BAM index not found' in err


# ---------------------------------------------------------------------------
# check subcommand — argument parsing
# ---------------------------------------------------------------------------

@pytest.fixture
def check_required_args(tmp_path):
    """Minimal argv tokens satisfying every required flag of the 'check'
    subcommand: just --vcfs and --fai (no workdir/models/annotations at
    all -- 'check' never touches any of those)."""
    fai = _touch(tmp_path / "ref.fa.fai")
    vcf = _touch(tmp_path / "manta.vcf")
    return {"fai": fai, "vcf": vcf}


@pytest.mark.unit
def test_check_subcommand_parses_fields(check_required_args):
    args = parsing.parse_arguments(
        arguments=[
            "check",
            "--vcfs", f"manta={check_required_args['vcf']}",
            "--fai", check_required_args["fai"],
            "--sample", "mysample",
            "--chrom", "chr1", "chr2",
            "--out", "/tmp/report.tsv",
        ]
    )
    assert args.command == "check"
    assert args.vcfs == [["manta", check_required_args["vcf"]]]
    assert args.fai == check_required_args["fai"]
    assert args.sample == "mysample"
    assert args.chrom == ["chr1", "chr2"]
    assert args.out == "/tmp/report.tsv"


@pytest.mark.unit
def test_check_subcommand_defaults(check_required_args):
    # Only --vcfs/--fai are required; --sample/--out default to None and
    # --chrom defaults to 'all' (same convention as 'call'/'multi').
    args = parsing.parse_arguments(
        arguments=[
            "check",
            "--vcfs", f"manta={check_required_args['vcf']}",
            "--fai", check_required_args["fai"],
        ]
    )
    assert args.sample is None
    assert args.chrom == "all"
    assert args.out is None


@pytest.mark.unit
@pytest.mark.parametrize("missing", ["--vcfs", "--fai"])
def test_check_missing_required_flag_exits(check_required_args, missing):
    full_args = [
        "check",
        "--vcfs", f"manta={check_required_args['vcf']}",
        "--fai", check_required_args["fai"],
    ]
    idx = full_args.index(missing)
    del full_args[idx:idx + 2]
    with pytest.raises(SystemExit):
        parsing.parse_arguments(arguments=full_args)


@pytest.mark.unit
def test_check_subcommand_has_no_annotation_or_model_flags(check_required_args):
    # 'check' reports on input VCFs only -- it has none of the annotation-
    # file, --models, --workdir, --bam or --pop flags 'call'/'multi' have.
    args = parsing.parse_arguments(
        arguments=[
            "check",
            "--vcfs", f"manta={check_required_args['vcf']}",
            "--fai", check_required_args["fai"],
        ]
    )
    for flag in ('annot_dir', 'models', 'workdir', 'bam', 'pop', 'pop_catalog',
                 'repeats', 'cgis', 'centromeres', 'gaps', 'althaps', 'vntrs',
                 'strs', 'gc'):
        assert not hasattr(args, flag)


@pytest.mark.unit
def test_check_never_touches_annotation_store(check_required_args, monkeypatch):
    # resolve_annotation_paths() early-returns for any command other than
    # 'call'/'multi', so parsing 'check' args must never call into the
    # managed-annotation-store machinery (which would otherwise prompt a
    # multi-GB first-use download).
    def _boom(*args, **kwargs):
        raise AssertionError('annotation store must not be touched for check')

    monkeypatch.setattr(parsing.annotations, 'annot_dir', _boom)
    monkeypatch.setattr(parsing.annotations, 'ensure_annotations', _boom)

    args = parsing.parse_arguments(
        arguments=[
            "check",
            "--vcfs", f"manta={check_required_args['vcf']}",
            "--fai", check_required_args["fai"],
        ]
    )
    assert args.command == "check"


# ---------------------------------------------------------------------------
# check subcommand — validate_inputs
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_validate_check_inputs_all_present_succeeds(check_required_args):
    args = parsing.parse_arguments(
        arguments=[
            "check",
            "--vcfs", f"manta={check_required_args['vcf']}",
            "--fai", check_required_args["fai"],
        ]
    )
    # Should not raise/exit, and (unlike 'call'/'multi') must not create a
    # workdir -- 'check' has none.
    parsing.validate_inputs(args)


@pytest.mark.unit
def test_validate_check_inputs_missing_fai_exits(check_required_args, capsys):
    args = parsing.parse_arguments(
        arguments=[
            "check",
            "--vcfs", f"manta={check_required_args['vcf']}",
            "--fai", check_required_args["fai"],
        ]
    )
    os.remove(args.fai)
    with pytest.raises(SystemExit):
        parsing.validate_inputs(args)
    err = capsys.readouterr().err
    assert '--fai file' in err
    assert args.fai in err


@pytest.mark.unit
def test_validate_check_inputs_missing_vcf_exits(check_required_args, capsys):
    args = parsing.parse_arguments(
        arguments=[
            "check",
            "--vcfs", f"manta={check_required_args['vcf']}",
            "--fai", check_required_args["fai"],
        ]
    )
    os.remove(check_required_args["vcf"])
    with pytest.raises(SystemExit):
        parsing.validate_inputs(args)
    err = capsys.readouterr().err
    assert 'VCF file for caller manta not found' in err


@pytest.mark.unit
def test_validate_check_inputs_no_vcfs_reports_required(check_required_args, capsys):
    args = parsing.parse_arguments(
        arguments=[
            "check",
            "--vcfs", f"manta={check_required_args['vcf']}",
            "--fai", check_required_args["fai"],
        ]
    )
    args.vcfs = []
    with pytest.raises(SystemExit):
        parsing.validate_inputs(args)
    err = capsys.readouterr().err
    assert '--vcfs is required' in err


@pytest.mark.unit
def test_validate_check_inputs_malformed_vcfs_token_exits(check_required_args, capsys):
    args = parsing.parse_arguments(
        arguments=[
            "check",
            "--vcfs", f"manta={check_required_args['vcf']}",
            "--fai", check_required_args["fai"],
        ]
    )
    args.vcfs = [["mantavcf"]]
    with pytest.raises(SystemExit):
        parsing.validate_inputs(args)
    err = capsys.readouterr().err
    assert 'Malformed --vcfs entry' in err
    assert 'mantavcf' in err


@pytest.mark.unit
def test_validate_check_inputs_rejects_duplicate_caller_labels(check_required_args, capsys):
    # Internal record ids are caller:ordinal, so two files under one label
    # would collide; every subcommand rejects that up front.
    args = parsing.parse_arguments(
        arguments=[
            "check",
            "--vcfs", f"delly={check_required_args['vcf']}", f"delly={check_required_args['vcf']}",
            "--fai", check_required_args["fai"],
        ]
    )
    with pytest.raises(SystemExit):
        parsing.validate_inputs(args)
    assert "Caller label used more than once: delly" in capsys.readouterr().err


@pytest.mark.unit
def test_validate_call_inputs_rejects_duplicate_caller_labels(call_required_args, capsys):
    args = parsing.parse_arguments(
        arguments=[
            "call", "--sample", "s",
            "--workdir", call_required_args["workdir"],
            "--fai", call_required_args["fai"],
            "--bam", call_required_args["bam"],
            "--annot-dir", call_required_args["annot_dir"],
            "--models", call_required_args["models_dir"],
            "--vcfs", f"delly={call_required_args['vcf']}", f"delly={call_required_args['vcf']}",
        ]
    )
    with pytest.raises(SystemExit):
        parsing.validate_inputs(args)
    assert "Caller label used more than once: delly" in capsys.readouterr().err


# ---------------------------------------------------------------------------
# build_parser — three subcommands
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_build_parser_exposes_all_three_subcommands():
    parser = parsing.build_parser()
    subparsers_action = next(
        action for action in parser._actions
        if isinstance(action, argparse._SubParsersAction)
    )
    assert set(subparsers_action.choices) == {'call', 'multi', 'check'}


@pytest.mark.unit
def test_parse_arguments_uses_build_parser(call_required_args):
    # parse_arguments() is a thin wrapper: build_parser() + parse_args() +
    # resolve_annotation_paths(). A namespace parsed via parse_arguments()
    # must be indistinguishable from one parsed by hand against build_parser().
    argv = [
        "call",
        "--sample", "mysample",
        "--workdir", call_required_args["workdir"],
        "--fai", call_required_args["fai"],
        "--bam", call_required_args["bam"],
        "--vcfs", f"manta={call_required_args['vcf']}",
    ]
    via_parse_arguments = parsing.parse_arguments(arguments=argv)
    via_build_parser = parsing.resolve_annotation_paths(parsing.build_parser().parse_args(argv))
    assert vars(via_parse_arguments) == vars(via_build_parser)
