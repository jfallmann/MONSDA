"""Long-read DTU slice: oarfish converter, DTU import modes, argument forwarding.

Portable by design: the R dependent checks skip when no usable ``Rscript`` is on
PATH, and the generator checks skip when snakemake is unavailable.  The R checks
exercise the real expressions taken from ``DRIMSEQ.R``/``DEXSEQ.R`` and the real
converter script, not re-implementations of them.
"""

import gzip
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "tests" / "data"
CONVERTER = REPO / "scripts" / "Analysis" / "oarfish_to_salmon.py"
DRIMSEQ = REPO / "scripts" / "Analysis" / "DTU" / "DRIMSEQ.R"
DEXSEQ = REPO / "scripts" / "Analysis" / "DTU" / "DEXSEQ.R"

RSCRIPT = shutil.which("Rscript")
needs_r = pytest.mark.skipif(
    RSCRIPT is None, reason="Rscript is not on PATH; activate an R environment"
)

# oarfish .quant with unequal transcript lengths but equal fractional counts
QUANT_ROWS = [("tx1", 1000, 10), ("tx2", 2000, 10), ("tx3", 500, 5)]


def _write_quant(path):
    lines = ["tname\tlen\tnum_reads"]
    lines += [f"{name}\t{length}\t{reads}" for name, length, reads in QUANT_ROWS]
    path.write_text("\n".join(lines) + "\n")
    return path


def _run_converter(inpath, outpath):
    return subprocess.run(
        [sys.executable, str(CONVERTER), "-i", str(inpath), "-o", str(outpath)],
        capture_output=True,
        text=True,
    )


def _read_quant_sf(path):
    opn = gzip.open if str(path).endswith(".gz") else open
    with opn(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    return header, rows


def _r(code, cwd=None):
    return subprocess.run(
        [RSCRIPT, "--no-save", "--no-restore", "-e", code],
        capture_output=True,
        text=True,
        cwd=cwd,
    )


def _r_has(pkg):
    if RSCRIPT is None:
        return False
    res = _r(f'cat(requireNamespace("{pkg}", quietly=TRUE))')
    return res.returncode == 0 and res.stdout.strip().endswith("TRUE")


def _call_at(text, start):
    """Return the balanced call starting at ``start`` in ``text``."""
    depth = 0
    for pos in range(start, len(text)):
        char = text[pos]
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
            if depth == 0:
                return text[start : pos + 1]
    raise AssertionError(f"unbalanced call at offset {start} in {text[start:start + 80]!r}")


def _import_expr(path, prefix, variable):
    """Extract the transcript-ID normalization call the script really uses."""
    text = path.read_text()
    start = text.index(prefix + variable + ")")
    call = _call_at(text, start)
    return call.replace(f"rownames({variable})", "ids")


# --------------------------------------------------------------------------- #
# oarfish -> salmon converter
# --------------------------------------------------------------------------- #


def test_converter_emits_salmon_columns(tmp_path):
    quant = _write_quant(tmp_path / "in.quant")
    out = tmp_path / "quant.sf"
    res = _run_converter(quant, out)
    assert res.returncode == 0, res.stderr

    header, rows = _read_quant_sf(out)
    assert header == ["Name", "Length", "EffectiveLength", "TPM", "NumReads"]
    assert [row[0] for row in rows] == [name for name, _, _ in QUANT_ROWS]
    # oarfish reports no effective length, so length is used for both columns
    for row, (_, length, reads) in zip(rows, QUANT_ROWS):
        assert float(row[1]) == length
        assert float(row[2]) == pytest.approx(length)
        assert float(row[4]) == pytest.approx(reads)
    assert sum(float(row[3]) for row in rows) == pytest.approx(1e6, rel=1e-6)


def test_converter_gzip_output(tmp_path):
    quant = _write_quant(tmp_path / "in.quant")
    out = tmp_path / "quant.sf.gz"
    res = _run_converter(quant, out)
    assert res.returncode == 0, res.stderr
    header, rows = _read_quant_sf(out)
    assert header[0] == "Name"
    assert len(rows) == len(QUANT_ROWS)


def test_converter_rejects_missing_column(tmp_path):
    quant = tmp_path / "bad.quant"
    quant.write_text("tname\tlen\n" + "tx1\t1000\n")
    res = _run_converter(quant, tmp_path / "out.sf")
    assert res.returncode != 0
    assert "num_reads" in res.stderr


# --------------------------------------------------------------------------- #
# import expressions taken from the real R scripts
# --------------------------------------------------------------------------- #


@needs_r
def test_drimseq_import_expression_strips_version_suffix():
    expr = _import_expr(DRIMSEQ, 'gsub("\\\\.[0-9]+$", "", rownames(', "cts_all")
    assert expr == 'gsub("\\\\.[0-9]+$", "", ids)'
    res = _r(
        'ids <- c("ENST00000456328.2", "ENST00000456328.2.1", "plain")\n'
        f"cat({expr}, sep='\\n')\n"
    )
    assert res.returncode == 0, res.stderr
    assert res.stdout.split() == ["ENST00000456328", "ENST00000456328.2", "plain"]


@needs_r
def test_dexseq_import_expression_strips_version_and_pipe_suffix():
    expr = _import_expr(DEXSEQ, 'sub("\\\\.[0-9]+$", "", sub("\\\\|.*", "", rownames(', "txi$counts")
    assert expr == 'sub("\\\\.[0-9]+$", "", sub("\\\\|.*", "", ids))'
    res = _r(
        'ids <- c("ENST00000456328.2|ENSG00000237613.2", "ENST00000456328.2", "plain")\n'
        f"cat({expr}, sep='\\n')\n"
    )
    assert res.returncode == 0, res.stderr
    assert res.stdout.split() == ["ENST00000456328", "ENST00000456328", "plain"]


def test_scripts_stop_on_ambiguous_normalized_ids():
    for script in (DRIMSEQ, DEXSEQ):
        text = script.read_text()
        assert "!duplicated(rownames(" not in text, script
        assert "ambiguous transcript IDs after version stripping" in text, script


# --------------------------------------------------------------------------- #
# import modes: exact ("no") versus legacy scaled ("scaledTPM")
# --------------------------------------------------------------------------- #


@needs_r
@pytest.mark.skipif(not _r_has("tximport"), reason="tximport is not installed")
def test_import_modes_exact_versus_scaled(tmp_path):
    quant = _write_quant(tmp_path / "in.quant")
    salmon_sf = tmp_path / "quant.sf"
    assert _run_converter(quant, salmon_sf).returncode == 0

    # a salmon style quant.sf with identical numbers, to show the converter
    # output is consumed exactly like salmon output
    salmon_like = tmp_path / "salmon_quant.sf"
    salmon_like.write_text(salmon_sf.read_text())

    samples = tmp_path / "samples"
    samples.mkdir()
    for name, src in (("conv", salmon_sf), ("salmon", salmon_like)):
        target = samples / name
        target.mkdir()
        shutil.copy(src, target / "quant.sf")

    code = f"""
suppressPackageStartupMessages(library(tximport))
files <- c(conv = "{samples}/conv/quant.sf", salmon = "{samples}/salmon/quant.sf")
for (mode in c("no", "scaledTPM")) {{
    txi <- tximport(files, type = "salmon", txOut = TRUE, countsFromAbundance = mode)
    m <- txi$counts
    for (tx in rownames(m)) {{
        cat("ROW", mode, tx, m[tx, "conv"], m[tx, "salmon"], "\\n")
    }}
}}
"""
    res = _r(code)
    assert res.returncode == 0, res.stderr

    counts = {}
    for line in res.stdout.splitlines():
        parts = line.split()
        if parts and parts[0] == "ROW":
            counts[(parts[1], parts[2])] = (float(parts[3]), float(parts[4]))

    # exact counts reproduce the reported NumReads
    assert counts[("no", "tx1")] == (10.0, 10.0)
    assert counts[("no", "tx2")] == (10.0, 10.0)
    assert counts[("no", "tx3")] == (5.0, 5.0)
    # scaled counts are redistributed by effective length: tx2 is twice as long
    # as tx1 at the same fractional count, so it loses half of its counts
    assert counts[("scaledTPM", "tx1")] == (10.0, 10.0)
    assert counts[("scaledTPM", "tx2")] == (5.0, 5.0)
    assert counts[("scaledTPM", "tx3")] == (10.0, 10.0)
    # converter output is consumed exactly like salmon output
    for mode in ("no", "scaledTPM"):
        for tx in ("tx1", "tx2", "tx3"):
            conv, salmon = counts[(mode, tx)]
            assert conv == salmon


# --------------------------------------------------------------------------- #
# invalid import mode
# --------------------------------------------------------------------------- #


@needs_r
@pytest.mark.parametrize("script", [DRIMSEQ, DEXSEQ], ids=["drimseq", "dexseq"])
def test_invalid_import_mode_stops(script, tmp_path):
    res = subprocess.run(
        [
            RSCRIPT,
            "--no-save",
            "--no-restore",
            str(script),
            "anno.gz",
            "anno.gtf",
            str(tmp_path),
            "T1-VS-T2",
            "none",
            "1",
            "min_samps_feature_expr = 1",
            "bogus",
        ],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert res.returncode != 0
    assert "unsupported DTU import mode 'bogus'" in res.stderr


# --------------------------------------------------------------------------- #
# generated argument forwarding (snakemake and nextflow)
# --------------------------------------------------------------------------- #

NF_EXPECTED = {
    "drimseq": (
        "params.gCOUNTTOOL = (get_always('drimseq_DTU_params_COUNTTOOL') ?: get_always('DTUCOUNTTOOL') ?: 'salmon').trim()",
        "params.gCOUNTPARAMS = get_always('drimseq_DTU_params_QUANT') ?: ''",
        "params.gDTUIMPORT = (params.gCOUNTTOOL == 'oarfish') ? "
        "((get_always('drimseq_DTU_params_DTU_IMPORT') ?: 'no').trim()) : 'scaledTPM'",
    ),
    "dexseq": (
        "params.gCOUNTTOOL = (get_always('dexseq_DTU_params_COUNTTOOL') ?: get_always('DTUCOUNTTOOL') ?: 'salmon').trim()",
        "params.gCOUNTPARAMS = get_always('dexseq_DTU_params_QUANT') ?: ''",
        "params.gDTUIMPORT = (params.gCOUNTTOOL == 'oarfish') ? "
        "((get_always('dexseq_DTU_params_DTU_IMPORT') ?: 'no').trim()) : 'scaledTPM'",
    ),
    "spit": (
        "params.gCOUNTTOOL = (get_always('spit_DTU_params_COUNTTOOL') ?: get_always('DTUCOUNTTOOL') ?: 'salmon').trim()",
        "params.gCOUNTPARAMS = get_always('spit_DTU_params_QUANT') ?: ''",
    ),
}

SMK_EXPECTED = {
    "drimseq": (
        "importmode = dtuimport if (COUNTTOOL == 'oarfish' and dtuimport) else "
        "('no' if COUNTTOOL == 'oarfish' else 'scaledTPM')",
        "importmode = importmode",
        "\\'{params.importmode}\\'",
    ),
    "dexseq": (
        "importmode = dtuimport if (COUNTTOOL == 'oarfish' and dtuimport) else "
        "('no' if COUNTTOOL == 'oarfish' else 'scaledTPM')",
        "importmode = importmode",
        "\\'{params.importmode}\\'",
    ),
}


@pytest.fixture
def dtu_workdir(tmp_path):
    """A temp working directory exposing tests/data inputs via symlinks."""
    for item in os.listdir(DATA):
        os.symlink(DATA / item, tmp_path / item)
    old = os.getcwd()
    os.chdir(tmp_path)
    try:
        yield tmp_path
    finally:
        os.chdir(old)


def _generate(engine, options, workdir):
    """Generate the DTU postflows and return (subdir, tool param strings)."""
    pytest.importorskip("snakemake")
    import MONSDA.Params as mp
    import MONSDA.Workflows as mw

    old_wf, old_env = mw.workflowpath, mw.envpath
    old_norm = mw.normalize_container_version
    mw.workflowpath = str(REPO / "workflows")
    mw.envpath = str(REPO / "envs") + os.sep
    mw.normalize_container_version = lambda version=None: "VERSION"
    try:
        with open(DATA / "config_DTU_test.json") as fh:
            cfg = json.load(fh)
        for tool, opts in options.items():
            cfg["DTU"]["Test"][tool]["OPTIONS"].update(opts)

        conditions = mp.get_conditions(cfg)
        _pre, sub, _post = mw.get_processes(cfg)
        combos = mw.get_combo(sub, cfg, conditions)
        samples = mp.get_samples_postprocess(cfg, "DTU")

        captured = []
        original = mw.nf_tool_params

        def spy(*a, **kw):
            out = original(*a, **kw)
            captured.append((a[4], out))
            return out

        mw.nf_tool_params = spy
        try:
            if engine == "nf":
                subdir = "SubPostFlows"
                mp.create_skeleton(subdir, None)
                mw.nf_make_post(
                    "DTU", cfg, samples, conditions, subdir, "INFO", combinations=combos
                )
            else:
                subdir = "SubPostSnakes"
                mp.create_skeleton(subdir, None)
                mw.make_post(
                    "DTU", cfg, samples, conditions, subdir, "INFO", combinations=combos
                )
        finally:
            mw.nf_tool_params = original
        return subdir, captured
    finally:
        mw.workflowpath, mw.envpath = old_wf, old_env
        mw.normalize_container_version = old_norm


def _subflow(workdir, subdir, tool, suffix):
    hits = sorted(Path(workdir, subdir).glob(f"*_DTU_{tool}_sub{suffix}"))
    assert len(hits) == 1, hits
    return hits[0].read_text()


@pytest.mark.parametrize("tool", ["drimseq", "dexseq", "spit"])
def test_nextflow_forwards_import_mode(tool, dtu_workdir):
    subdir, captured = _generate(
        "nf", {tool: {"COUNTTOOL": "oarfish", "DTU_IMPORT": "scaledTPM"}}, dtu_workdir
    )
    text = _subflow(dtu_workdir, subdir, tool, "flow.nf")
    for line in NF_EXPECTED[tool]:
        assert line in text, line

    params = dict(captured)
    assert f"--{tool}_DTU_params_COUNTTOOL 'oarfish '" in params[f"{tool}_DTU"]

    if tool == "spit":
        # SPIT.py takes no import mode argument
        assert "gDTUIMPORT" not in text
        assert "$imp" not in text
    else:
        assert "imp = \"'${params.gDTUIMPORT}'\"" in text
        assert "$dparams $imp &> log &&" in text


@pytest.mark.parametrize("tool", ["drimseq", "dexseq"])
def test_snakemake_forwards_import_mode(tool, dtu_workdir):
    subdir, _captured = _generate(
        "smk", {tool: {"COUNTTOOL": "oarfish", "DTU_IMPORT": "scaledTPM"}}, dtu_workdir
    )
    text = _subflow(dtu_workdir, subdir, tool, "snake.smk")
    for line in SMK_EXPECTED[tool]:
        assert line in text, line


@pytest.mark.parametrize("engine", ["nf", "smk"])
def test_import_mode_defaults_to_salmon_scaled_tpm(engine, dtu_workdir):
    subdir, captured = _generate(engine, {}, dtu_workdir)
    if engine == "nf":
        text = _subflow(dtu_workdir, subdir, "drimseq", "flow.nf")
        assert "get_always('drimseq_DTU_params_COUNTTOOL') ?: get_always('DTUCOUNTTOOL') ?: 'salmon'" in text
        assert "gDTUIMPORT" in text
        assert "scaledTPM" in text
        params = dict(captured)
        assert "--drimseq_DTU_params_COUNTTOOL" not in params["drimseq_DTU"]
    else:
        text = _subflow(dtu_workdir, subdir, "drimseq", "snake.smk")
        assert "COUNTTOOL', 'salmon'" in text
        assert "importmode = importmode" in text


@pytest.mark.parametrize("engine", ["nf", "smk"])
def test_import_mode_ignored_for_salmon(engine, dtu_workdir):
    """DTU_IMPORT without COUNTTOOL=oarfish must stay on the short read default."""
    subdir, _captured = _generate(
        "nf" if engine == "nf" else "smk", {"drimseq": {"DTU_IMPORT": "no"}}, dtu_workdir
    )
    if engine == "nf":
        text = _subflow(dtu_workdir, subdir, "drimseq", "flow.nf")
        assert "((get_always('drimseq_DTU_params_DTU_IMPORT') ?: 'no').trim()) : 'scaledTPM'" in text
    else:
        text = _subflow(dtu_workdir, subdir, "drimseq", "snake.smk")
        assert "('no' if COUNTTOOL == 'oarfish' else 'scaledTPM')" in text


@pytest.mark.parametrize(
    "options,expected,used",
    [
        ({"QUANT": "--numReads 5"}, {"QUANT": "--numReads 5"}, None),
        ({"COUNT": "--numReads 5"}, {"QUANT": "--numReads 5"}, "OPTIONS.COUNT"),
        (
            {"QUANT": "--numReads 9", "COUNT": "--numReads 5"},
            {"QUANT": "--numReads 9"},
            None,
        ),
        (
            {"QUANT": "--numReads 9", "COUNT": "  "},
            {"QUANT": "--numReads 9"},
            "OPTIONS.QUANT",
        ),
        (
            {"QUANT": "", "COUNT": "--numReads 5"},
            {"QUANT": "--numReads 5"},
            "OPTIONS.COUNT",
        ),
        ({"QUANT": " ", "COUNT": ""}, {"QUANT": " "}, None),
        ({"QUANT": " --numReads 5 "}, {"QUANT": " --numReads 5 "}, None),
        ({"INDEX": "-d GENOMES/idx"}, {"INDEX": "-d GENOMES/idx"}, None),
    ],
    ids=[
        "quant-only",
        "count-only",
        "both-set",
        "count-blank",
        "quant-blank",
        "both-blank",
        "verbatim",
        "unrelated-only",
    ],
)
def test_count_option_warning(options, expected, used, caplog):
    import MONSDA.Params as mp

    original = dict(options)
    with caplog.at_level("WARNING"):
        result = mp.normalize_count_options(options)
    assert options == original
    assert dict(result) == expected
    warnings = [r.message for r in caplog.records if r.levelname == "WARNING"]
    if used is None:
        assert warnings == []
    else:
        assert len(warnings) == 1
        assert "OPTIONS.QUANT" in warnings[0]
        assert "OPTIONS.COUNT" in warnings[0]
        assert warnings[0].endswith(used)
        assert "--numReads" not in warnings[0]


COUNT_STAGE_TOOLS = [
    ("COUNTING", "countreads", {}),
    ("DE", "deseq2", {}),
    ("DE", "deseq2", {"COUNTTOOL": "oarfish"}),
    ("DEU", "edger", {}),
    ("DAS", "edger", {}),
    ("DTU", "drimseq", {}),
    ("DTU", "drimseq", {"COUNTTOOL": "oarfish"}),
]

COUNT_OPTION_CASES = [
    ({"QUANT": "--numReads 5"}, "--numReads 5"),
    ({"COUNT": "--numReads 5"}, "--numReads 5"),
    ({"QUANT": "--numReads 9", "COUNT": "--numReads 5"}, "--numReads 9"),
    ({"QUANT": "  ", "COUNT": "--numReads 5"}, "--numReads 5"),
    ({"QUANT": "--numReads 9", "COUNT": ""}, "--numReads 9"),
    ({"QUANT": " --numReads 5 "}, " --numReads 5 "),
]


@pytest.mark.parametrize(
    "stage,tool,extra",
    COUNT_STAGE_TOOLS,
    ids=["-".join([s, t] + sorted(e.values())) for s, t, e in COUNT_STAGE_TOOLS],
)
@pytest.mark.parametrize(
    "options,expected",
    COUNT_OPTION_CASES,
    ids=["quant", "count", "both", "count-blank", "quant-blank", "verbatim"],
)
def test_count_options_reach_both_engines(
    stage, tool, extra, options, expected, dtu_workdir
):
    pytest.importorskip("snakemake")
    import MONSDA.Params as mp
    import MONSDA.Workflows as mw

    with open(DATA / "config_Test.json") as fh:
        cfg = json.load(fh)
    cfg[stage]["Test"][tool]["OPTIONS"] = dict(options, **extra)
    sample = mp.get_samples_postprocess(cfg, stage)[0]

    opts = dict(mp.tool_params(sample, None, cfg, stage, tool)["OPTIONS"])
    assert opts["QUANT"] == expected
    assert "COUNT" not in opts

    nf = mw.nf_tool_params(sample, None, cfg, stage, tool, tool, None, None)
    assert f"--{tool}_params_QUANT '{expected} '" in nf
    assert "_params_COUNT '" not in nf
