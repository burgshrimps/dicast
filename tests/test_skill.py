"""Keep the agent skill (.agents/skills/dicast, mirrored into .claude/skills)
in sync with the CLI: the frontmatter must parse and every --flag the skill
mentions must exist in the argparse parser, so the skill cannot silently go
stale when a flag is renamed or removed.

Adapted from cuban's tests/test_skill.py. dicast's parser is split into a
top-level parser plus three subparsers ('call', 'multi', 'check'), each with
its own flags, so -- unlike cuban's single flat parser -- flags are collected
from build_parser() *and* every subparser it defines.
"""

import argparse
import re
from pathlib import Path

from dicast.parsing import build_parser

REPO_ROOT = Path(__file__).resolve().parent.parent
SKILL_DIR = REPO_ROOT / ".agents" / "skills" / "dicast"
SKILL_MD = SKILL_DIR / "SKILL.md"
CLAUDE_LINK = REPO_ROOT / ".claude" / "skills" / "dicast"

# Tools the skill's "Fixing inputs" recipes shell out to; their flags are not
# dicast's and should not be checked against dicast's parser.
EXTERNAL_TOOLS = ("bcftools", "samtools", "tabix", "Rscript", "gridss", "cnvnator2VCF", "cuban")


def _frontmatter(text):
    match = re.match(r"^---\n(.*?)\n---\n", text, re.S)
    assert match, "SKILL.md must start with a YAML frontmatter block"
    fields = {}
    for line in match.group(1).splitlines():
        key, sep, value = line.partition(":")
        assert sep, f"malformed frontmatter line: {line!r}"
        fields[key.strip()] = value.strip()
    return fields


def _all_parsers():
    """Yields the top-level parser plus every subparser it defines."""
    parser = build_parser()
    yield parser
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            yield from action.choices.values()


def _all_flags():
    """Every --flag/-x option string defined anywhere in dicast's parser
    (top level and every subcommand)."""
    flags = set()
    for parser in _all_parsers():
        for action in parser._actions:
            flags.update(action.option_strings)
    return flags


def test_frontmatter_has_name_and_description():
    fields = _frontmatter(SKILL_MD.read_text())
    assert fields["name"] == "dicast"
    # Both Claude Code and Codex cap the description at 1024 characters.
    assert 0 < len(fields["description"]) <= 1024


def test_every_flag_in_skill_exists_in_cli():
    parser_flags = _all_flags()
    text = "\n".join(p.read_text() for p in SKILL_DIR.rglob("*.md"))
    # Flags in the "Fixing inputs" shell-out recipes belong to external
    # tools, not dicast; ignore lines mentioning them.
    lines = [
        line for line in text.splitlines()
        if not any(tool in line for tool in EXTERNAL_TOOLS)
    ]
    # Allows underscores (dicast has --sv_types, --pop_catalog-style flags).
    mentioned = set(re.findall(r"(?<![\w-])(--[a-z][a-z_-]*)", "\n".join(lines)))
    unknown = sorted(mentioned - parser_flags)
    assert not unknown, f"SKILL.md mentions flags the CLI does not have: {unknown}"


def test_claude_skill_mirrors_agents_skill():
    assert CLAUDE_LINK.is_symlink(), ".claude/skills/dicast must be a symlink"
    assert CLAUDE_LINK.resolve() == SKILL_DIR.resolve()
    assert (CLAUDE_LINK / "SKILL.md").is_file()
