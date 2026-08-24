"""Unit tests for ``dicast.annotations`` — the managed annotation store.

Network is never touched: download behavior is exercised only up to the
decision to fetch (``missing_files`` / ``ensure_annotations`` no-op paths).
"""
import pytest

from dicast import annotations


@pytest.mark.unit
def test_env_var_overrides_annot_dir(monkeypatch, tmp_path):
    monkeypatch.setenv("DICAST_DATA_DIR", str(tmp_path))
    assert annotations.annot_dir() == tmp_path


@pytest.mark.unit
def test_stored_config_choice_wins_without_env(monkeypatch, tmp_path):
    monkeypatch.delenv("DICAST_DATA_DIR", raising=False)
    chosen = tmp_path / "chosen_store"
    monkeypatch.setattr(annotations, "_load_config", lambda: {"annot_dir": str(chosen)})
    assert annotations.annot_dir() == chosen


@pytest.mark.unit
def test_suggested_dir_is_repo_annot_for_checkout(monkeypatch, tmp_path):
    monkeypatch.delenv("DICAST_DATA_DIR", raising=False)
    monkeypatch.setattr(annotations, "_load_config", lambda: {})
    # Running from this repo checkout, the suggestion is <repo>/annot.
    assert annotations._suggested_annot_dir().name == "annot"


@pytest.mark.unit
def test_missing_files_reports_only_absent(tmp_path):
    present = "hg38_centromeres.tsv"
    (tmp_path / present).write_text("")
    missing = annotations.missing_files(tmp_path)
    assert present not in missing
    assert set(missing) == set(annotations.ANNOT_FILES) - {present}


@pytest.mark.unit
def test_ensure_annotations_noop_when_populated(tmp_path):
    for name in annotations.ANNOT_FILES:
        (tmp_path / name).write_text("")
    # Fully populated store: returns without attempting any download.
    annotations.ensure_annotations(list(annotations.ANNOT_FILES), tmp_path)


@pytest.mark.unit
def test_fetch_checksum_mismatch_rejects_file(monkeypatch, tmp_path):
    # A downloaded file whose checksum does not match must be discarded.
    import io
    import urllib.request

    monkeypatch.setattr(urllib.request, "urlopen", lambda url: io.BytesIO(b"corrupt"))
    with pytest.raises(SystemExit, match="checksum"):
        annotations.fetch_annotations(["hg38_centromeres.tsv"], dest=tmp_path, quiet=True)
    assert not (tmp_path / "hg38_centromeres.tsv").exists()
    assert not (tmp_path / "hg38_centromeres.tsv.part").exists()
