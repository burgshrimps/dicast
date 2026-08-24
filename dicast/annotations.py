"""Locating, storing and fetching the hg38 reference annotations."""

import hashlib
import json
import os
import shutil
import sys
import urllib.request
from pathlib import Path

RELEASE_URL_BASE = "https://github.com/burgshrimps/dicast/releases/latest/download/"

# filename -> (md5, approximate size) of every annotation asset on the release.
ANNOT_FILES = {
    "hg38_repeatmasker.tsv": ("1210bf55914dc3f5ac6de7fe471efffb", "460 MB"),
    "hg38_cpg_islands.tsv": ("138b61a7ea544d827018261828ff6672", "2 MB"),
    "hg38_centromeres.tsv": ("f87165e17ab6b67c1f31c35326b34401", "4 KB"),
    "hg38_asmb_gaps.tsv": ("3d71d0d2e35fa1dabd7ffee29950e142", "44 KB"),
    "hg38_alt_haps.tsv": ("7fea41b4dd1deb91821ddef9e314ee66", "1.5 MB"),
    "hg38_vntrs_chaisson.bed": ("688d7c0400fcb3ff14a865ac24ef14d3", "830 KB"),
    "hg38_strs_chaisson.bed": ("a922016dde8a814e71482beb25afb97f", "27 MB"),
    "hg38_gc_content.bw": ("40fd8cc989e45eeab0d992aad255e3aa", "1.7 GB"),
}

TOTAL_SIZE = "~2.1 GB"


def _repo_root():
    """The repo checkout the package runs from, or None for a site-packages install."""
    root = Path(__file__).resolve().parent.parent
    if (root / "pyproject.toml").is_file() and os.access(root, os.W_OK):
        return root
    return None


CONFIG_PATH = Path.home() / ".dicast" / "config.json"


def _load_config():
    try:
        with open(CONFIG_PATH) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return {}


def _save_config(config):
    CONFIG_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(CONFIG_PATH, "w") as fh:
        json.dump(config, fh, indent=2)


def _suggested_annot_dir():
    """The prefilled storage suggestion: annot/ in the repo root when running
    from a checkout, ~/.dicast otherwise."""
    root = _repo_root()
    if root is not None:
        return root / "annot"
    return Path.home() / ".dicast"


def annot_dir(ask=False):
    """Where the annotation files are stored: $DICAST_DATA_DIR if set, else the
    location the user chose on first use (persisted in ~/.dicast/config.json).

    With ask=True and no stored choice yet, the user is prompted once for the
    location (prefilled with the suggestion; Enter confirms). In
    non-interactive sessions the suggestion is used without asking.
    """
    env_dir = os.environ.get("DICAST_DATA_DIR")
    if env_dir:
        return Path(env_dir)
    config = _load_config()
    stored = config.get("annot_dir")
    if stored:
        return Path(stored)
    suggestion = _suggested_annot_dir()
    if not ask:
        return suggestion
    if sys.stdin.isatty() and sys.stderr.isatty():
        print(
            f"dicast stores the hg38 reference annotations ({TOTAL_SIZE}, "
            "downloaded once) that its feature engineering depends on.",
            file=sys.stderr,
        )
        answer = input(f"Storage location [{suggestion}]: ").strip()
        chosen = Path(answer).expanduser() if answer else suggestion
    else:
        chosen = suggestion
        print(f"[dicast] storing the annotation files in {chosen} "
              "(set DICAST_DATA_DIR to change)", file=sys.stderr)
    try:
        chosen.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        raise SystemExit(f"cannot create annotation directory {chosen}: {e}") from e
    config["annot_dir"] = str(chosen)
    _save_config(config)
    return chosen


def _md5(path):
    digest = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def missing_files(dest, filenames=None):
    """The subset of filenames (default: all) not present in dest."""
    filenames = filenames or list(ANNOT_FILES)
    dest = Path(dest)
    return [name for name in filenames if not (dest / name).is_file()]


def fetch_annotations(filenames=None, dest=None, quiet=False):
    """Download the given annotation files (default: all) to dest (default:
    the annot directory, asking once where that should be)."""
    filenames = filenames or list(ANNOT_FILES)
    dest = Path(dest) if dest else annot_dir(ask=True)
    dest.mkdir(parents=True, exist_ok=True)
    for name in filenames:
        md5_expected, size = ANNOT_FILES[name]
        target = dest / name
        if target.is_file():
            continue
        url = RELEASE_URL_BASE + name
        tmp = target.with_suffix(target.suffix + ".part")
        if not quiet:
            print(f"downloading {url} ({size})\n        -> {target}", file=sys.stderr)
        try:
            with urllib.request.urlopen(url) as resp, open(tmp, "wb") as out:
                shutil.copyfileobj(resp, out)
        except OSError as e:
            tmp.unlink(missing_ok=True)
            raise SystemExit(f"download of {name} failed: {e}") from e
        md5_actual = _md5(tmp)
        if md5_actual != md5_expected:
            tmp.unlink(missing_ok=True)
            raise SystemExit(
                f"downloaded {name} failed its checksum "
                f"(expected {md5_expected}, got {md5_actual}); "
                f"download it manually from {url}")
        tmp.replace(target)
        if not quiet:
            print(f"done ({target.stat().st_size / 1e6:.0f} MB)", file=sys.stderr)
    return dest


def ensure_annotations(filenames, dest):
    """Fetch any of filenames missing from dest, announcing the one-time download."""
    needed = missing_files(dest, filenames)
    if not needed:
        return
    print(
        f"[dicast] {len(needed)} annotation file(s) not found in {dest}; "
        "downloading them (one time)...", file=sys.stderr,
    )
    try:
        fetch_annotations(needed, dest=dest, quiet=False)
    except SystemExit as e:
        raise SystemExit(
            f"{e}\nRun 'dicast-fetch-annotations' once you are online, or "
            "point --annot-dir (or the individual annotation flags) at your "
            "own copies.") from e


def fetch_main(argv=None):
    import argparse

    p = argparse.ArgumentParser(
        prog="dicast-fetch-annotations",
        description="Download the hg38 reference annotations used by dicast "
                    f"({TOTAL_SIZE}) to {annot_dir()}/ (override with --dest "
                    "or the DICAST_DATA_DIR environment variable).",
    )
    p.add_argument("--dest", help="write the files to this directory instead")
    args = p.parse_args(argv)
    fetch_annotations(dest=args.dest)
