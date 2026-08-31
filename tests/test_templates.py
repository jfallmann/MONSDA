"""Static convention checks on the workflow templates.

These do not need snakemake, nextflow or any tool binary, they only read the
templates in ``workflows/``. Each check corresponds to a class of bug that has
already shipped at least once, see the docstrings of the individual tests.
"""

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
WORKFLOWS = REPO / "workflows"

UPSTREAM = "(?:MAPPED|TRIMMED_FASTQ|DEDUP_FASTQ|UNMAPPED)"
POSTDIR = re.compile(r'"(?:COUNTS|DE|DEU|DAS|DTU|PEAKS|TRACKS|FUSIONS|CIRCS|BED)/')

SECTION = re.compile(r"^\s*(\w+):")
SECTION_KEYS = {
    "input", "output", "log", "params", "conda", "container", "threads",
    "shell", "run", "script", "message", "priority", "resources", "benchmark",
    "group", "retries", "wildcard_constraints", "shadow", "envmodules",
    "cache", "default_target", "localrule", "wrapper", "notebook",
}
EXPAND = re.compile(r'expand\(\s*(?:f?")([^"]*)"\s*,\s*([^)]*)\)')
RULEHEAD = re.compile(r"^(\s*)(?:rule|checkpoint)\s+(\w+)\s*:\s*$")

NFWORKFLOW = re.compile(
    r"^workflow\s+(?:COUNTING|DE|DEU|DAS|DTU|PEAKS|TRACKS|FUSIONS|CIRCS)\s*\{"
)
NFUPSTREAM = re.compile(UPSTREAM + r"/\$\{SCOMBO\}")


def _postprocess_smk():
    """Templates that write into a postprocessing result directory."""
    return sorted(p for p in WORKFLOWS.glob("*.smk") if POSTDIR.search(p.read_text()))


def _postprocess_nf():
    """Templates that define a postprocessing sub-workflow."""
    out = []
    for path in sorted(WORKFLOWS.glob("*.nf")):
        lines = path.read_text().split("\n")
        start = next((i for i, l in enumerate(lines) if NFWORKFLOW.match(l)), None)
        if start is not None:
            out.append((path, lines, start))
    return out


def test_postprocess_smk_inputs_use_scombo():
    """Upstream inputs of a postprocessing rule live under the mapping combo.

    In a postprocessing sub-workflow ``combo`` is ``<mappingcombo>/<tool>`` and
    ``scombo`` is ``<mappingcombo>``, so anything read back from ``MAPPED/``,
    ``TRIMMED_FASTQ/`` or ``DEDUP_FASTQ/`` has to be expanded with ``scombo``.
    Using ``{combo}`` there looks for the files under a subdirectory named after
    the postprocessing tool, which never exists.
    """
    bad = []
    for path in _postprocess_smk():
        section = None
        for lineno, line in enumerate(path.read_text().split("\n"), 1):
            if line.lstrip().startswith("#"):
                continue
            match = SECTION.match(line)
            if match and match.group(1) in SECTION_KEYS:
                section = match.group(1)
            elif line and not line[0].isspace():
                section = None
            if section != "input":
                continue
            for call in EXPAND.finditer(line):
                pattern, kwargs = call.group(1), call.group(2)
                for placeholder in re.findall(UPSTREAM + r"/\{(s?combo)\}", pattern):
                    bound = re.search(r"\b%s\s*=\s*(\w+)" % placeholder, kwargs)
                    if (bound.group(1) if bound else None) != "scombo":
                        bad.append(
                            f"{path.name}:{lineno}: upstream {{{placeholder}}} is not "
                            f"bound to scombo"
                        )
            if re.search(UPSTREAM + r"/\{combo\}", EXPAND.sub("", line)):
                bad.append(
                    f"{path.name}:{lineno}: upstream path uses {{combo}}, expand it "
                    f"with scombo instead: {line.strip()[:80]}"
                )
    assert not bad, "postprocessing inputs read from the wrong combo:\n" + "\n".join(bad)


def test_smk_rules_have_a_body():
    """A rule header directly followed by another header is a copy/paste slip.

    ``annotatebed.smk`` carried a duplicated ``rule BamToBed:`` line, which the
    snakemake parser happily accepts, silently dropping one of the two.
    """
    bad = []
    for path in sorted(WORKFLOWS.glob("*.smk")):
        previous = None
        for lineno, line in enumerate(path.read_text().split("\n"), 1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            match = RULEHEAD.match(line)
            if match and previous is not None:
                bad.append(
                    f"{path.name}:{lineno}: rule {match.group(2)!r} directly follows "
                    f"rule header {previous!r}, the first one has no body"
                )
            previous = match.group(2) if match else None
    assert not bad, "rule without a body:\n" + "\n".join(bad)


def test_postprocess_nf_builds_its_own_channels():
    """Postprocessing nextflow workflows must not consume ``collection``.

    ``RunMONSDA`` calls every postprocessing sub-workflow as ``SUBWORK(dummy())``
    where ``dummy`` is just ``LOGS/MONSDA.log``. A workflow that filters or maps
    ``collection`` therefore ends up with an empty channel and silently runs on
    no sample at all. The upstream files have to be listed from ``*SAMPLES``.
    """
    bad = []
    for path, lines, start in _postprocess_nf():
        body = [l for l in lines[start:] if not re.match(r"^\s*take:\s*collection", l)]
        source = "\n".join(body)
        if re.search(r"\bcollection\b", source):
            bad.append(
                f"{path.name}: consumes 'collection', which only ever carries the "
                f"MONSDA log dummy"
            )
        if not re.search(r"\b\w*SAMPLES\b", source):
            bad.append(f"{path.name}: builds no input channel from a *SAMPLES list")
    assert not bad, "postprocessing workflow never reaches a sample:\n" + "\n".join(bad)


def test_postprocess_nf_upstream_paths_use_combo():
    """In nextflow the combo variables carry the opposite meaning.

    ``--COMBO`` is the mapping combo and ``--SCOMBO`` is
    ``<mappingcombo>/<tool>``, so upstream files are read from ``${COMBO}`` and
    only results are published under ``${SCOMBO}``.
    """
    bad = []
    for path, lines, start in _postprocess_nf():
        for lineno, line in enumerate(lines[start:], start + 1):
            if NFUPSTREAM.search(line):
                bad.append(
                    f"{path.name}:{lineno}: upstream path uses ${{SCOMBO}}, must be "
                    f"${{COMBO}}: {line.strip()[:80]}"
                )
    assert not bad, "upstream nextflow path uses the wrong combo:\n" + "\n".join(bad)


def test_checks_cover_all_templates():
    """Guard against the selectors silently matching nothing."""
    assert len(_postprocess_smk()) > 20
    assert len(_postprocess_nf()) > 20
