import glob
import gzip
import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest
import yaml

pytest.importorskip("snakemake")

from snakemake.common.configfile import load_configfile  # noqa: E402

import MONSDA.Params as mp  # noqa: E402
import MONSDA.Workflows as mw  # noqa: E402
from MONSDA.Containers import normalize_container_version  # noqa: E402
from scripts.container_matrix import SKIP_ENVIRONMENTS, environment_names  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "tests" / "data"
GOLDEN = REPO / "tests" / "golden"

UPDATE = bool(os.environ.get("UPDATE_GOLDEN"))


@pytest.fixture(autouse=True)
def _repo_template_paths(monkeypatch):
    """Point the generators at the in-repo workflow/env templates.

    In a dev checkout MONSDA is not pip-installed, so the module-level
    ``workflowpath``/``envpath`` (which assume a site-packages layout) do not
    resolve. Redirect them to the repository copies for the duration of a test.
    """
    old_wf, old_env = mw.workflowpath, mw.envpath
    mw.workflowpath = str(REPO / "workflows")
    mw.envpath = str(REPO / "envs") + os.sep
    monkeypatch.setattr(mw, "normalize_container_version", lambda version=None: "VERSION")
    yield
    mw.workflowpath, mw.envpath = old_wf, old_env


def _load(name):
    return load_configfile(str(DATA / name))


def _normalize(text, workdir):
    return (
        text.replace(str(workdir), "<WORKDIR>")
        .replace(mw.binpath, "<REPO>/MONSDA/MONSDA/scripts")
        .replace(str(REPO), "<REPO>")
    )


def _check_golden(relpath, content):
    path = GOLDEN / relpath
    if UPDATE:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
        return
    assert path.exists(), (
        f"missing golden snapshot {relpath!s}; regenerate with UPDATE_GOLDEN=1"
    )
    assert content == path.read_text(), (
        f"generator output drifted from golden snapshot {relpath!s}; "
        f"review the diff and, if intended, regenerate with UPDATE_GOLDEN=1"
    )


@pytest.fixture
def workdir(tmp_path):
    """A temp working directory exposing tests/data inputs via symlinks."""
    for item in os.listdir(DATA):
        os.symlink(DATA / item, tmp_path / item)
    old = os.getcwd()
    os.chdir(tmp_path)
    try:
        yield tmp_path
    finally:
        os.chdir(old)


# --------------------------------------------------------------------------- #
# Unit tests: pure helpers
# --------------------------------------------------------------------------- #


class TestCreateSubworkflow:
    def test_returns_tool_and_config(self):
        cfg = _load("config_fgumi_test.json")
        conditions = mp.get_conditions(cfg)
        toollist, configs = mw.create_subworkflow(cfg, "DEDUP", conditions)
        assert toollist == [["fgumi", "fgumi"]]
        assert len(configs) == 1
        assert configs[0]["RUNDEDUP"] == "enabled"
        assert "fgumi" in configs[0]["DEDUP"]["Test"]["umi"]

    def test_tools_key_without_env_bin(self):
        cfg = _load("config_Test.json")
        conditions = mp.get_conditions(cfg)
        toollist, configs = mw.create_subworkflow(cfg, "MAPPING", conditions)
        tools = sorted(t[0] for t in toollist)
        assert tools
        assert all(len(pair) == 2 for pair in toollist)

    def test_missing_env_bin_logs_error(self, caplog):
        cfg = {
            "MAXTHREADS": "1",
            "WORKFLOWS": "MAPPING",
            "SETTINGS": {"Test": {"umi": {"SAMPLES": ["S1"]}}},
            "MAPPING": {"Test": {"umi": {}}},
        }
        with caplog.at_level("ERROR"):
            toollist, _ = mw.create_subworkflow(cfg, "MAPPING", [("Test", "umi")])
        assert toollist == []
        assert any(
            "No tool environment found for MAPPING" in r.message
            for r in caplog.records
        )
        assert any(
            "No tool binary found for MAPPING" in r.message for r in caplog.records
        )


class TestGetCombo:
    def test_single_tool(self):
        cfg = _load("config_fgumi_test.json")
        conditions = mp.get_conditions(cfg)
        combos = mw.get_combo(["DEDUP"], cfg, conditions)
        materialized = {c: [list(x) for x in combos[c]] for c in combos}
        assert materialized == {("Test", "umi"): [[{"DEDUP": "fgumi"}]]}

    def test_none_when_no_workflows(self):
        cfg = _load("config_fgumi_test.json")
        conditions = mp.get_conditions(cfg)
        assert mw.get_combo(None, cfg, conditions) is None
        assert mw.get_combo([], cfg, conditions) is None

    def test_multi_tool_product(self):
        cfg = _load("config_Test.json")
        conditions = mp.get_conditions(cfg)
        combos = mw.get_combo(["MAPPING"], cfg, conditions)
        counts = {c: len(list(combos[c])) for c in combos}
        assert all(n >= 1 for n in counts.values())


class TestNfRwGetAlways:
    """Literal get_always('KEY')/get_always(\"KEY\") lookups must survive
    _nf_rw byte-exact; only dynamic get_always(KEY) and GString arguments
    are rewritten."""

    def test_plain_literal_protected(self):
        text = "PCOMBO = get_always('COMBO') ?: 'none'"
        assert mw._nf_rw(text, {"COMBO"}, set()) == text

    def test_double_quoted_literal_protected(self):
        text = 'PCOMBO = get_always("COMBO") ?: "none"'
        assert mw._nf_rw(text, {"COMBO"}, set()) == text

    def test_whitespace_and_escaping_preserved(self):
        text = "X = get_always( 'KEY\\'S' )"
        assert mw._nf_rw(text, {"KEY"}, set()) == text

    def test_surrounding_bare_and_interpolation_rewritten(self):
        text = "X = COMBO + get_always('COMBO') + $COMBO"
        assert mw._nf_rw(text, {"COMBO"}, set()) == (
            "X = params.gCOMBO + get_always('COMBO') + ${params.gCOMBO}"
        )

    def test_dynamic_args_rewritten(self):
        assert mw._nf_rw("Y = get_always(COMBO)", {"COMBO"}, set()) == (
            "Y = get_always(params.gCOMBO)"
        )

    def test_gstring_arg_rewritten(self):
        assert mw._nf_rw('Z = get_always("$COMBO")', {"COMBO"}, set()) == (
            'Z = get_always("${params.gCOMBO}")'
        )

    def test_template_literal_keys_preserved(self):
        """Every literal get_always lookup in a workflow template must stay
        byte-exact even when its key is itself a rewritten global name (the
        key/global intersection that broke PCOMBO)."""
        calls = set()
        for path in (REPO / "workflows").glob("*.nf"):
            for m in mw._GA_LITERAL.finditer(path.read_text()):
                if m.group(1) == '"' and "$" in m.group(2):
                    continue
                calls.add((path.name, m.group(0), m.group(2)))
        assert calls
        for name, call, key in calls:
            assert mw._nf_rw(call, {key}, set()) == call, (
                f"{name}: literal lookup {call!r} rewritten"
            )


def test_nf_make_workflow_pcombo_preserves_key(workdir):
    """The full header+template generation must keep the literal
    get_always('COMBO') lookup in PCOMBO byte-exact instead of rewriting
    the key to params.gCOMBO."""
    cfg = _load("config_DTU_test.json")
    conditions = mp.get_conditions(cfg)
    _pre, sub, _post = mw.get_processes(cfg)
    combos = mw.get_combo(sub, cfg, conditions)
    samples = mp.get_samples_postprocess(cfg, "DTU")
    subdir = "SubPostFlows"
    mp.create_skeleton(subdir, None)
    mw.nf_make_post(
        "DTU", cfg, samples, conditions, subdir, "INFO", combinations=combos
    )
    produced = sorted(
        os.path.relpath(p, workdir)
        for p in glob.glob(os.path.join(subdir, "*_subflow.nf"))
    )
    assert produced, "generator produced no subflows"
    for rel in produced:
        text = Path(workdir / rel).read_text()
        assert "params.gPCOMBO = get_always('COMBO') ?: 'none'" in text, rel
        assert "get_always('params.gCOMBO')" not in text, rel


def test_diego_nf_log_path_uses_declared_comparisons(workdir):
    cfg = _load("config_Test.json")
    cfg["DAS"]["TOOLS"] = {"diego": "Analysis/DAS/DIEGO.py"}
    cfg["DAS"]["Test"] = {
        "diego": {
            "OPTIONS": {
                "COUNT": "-t exon -g gene_id -f",
                "DAS": "",
            }
        }
    }
    cfg["DAS"]["COMPARABLE"] = {"T1-VS-T2": [["T1"], ["T2"]]}
    conditions = mp.get_conditions(cfg)
    _pre, sub, _post = mw.get_processes(cfg)
    combos = mw.get_combo(sub, cfg, conditions)
    samples = mp.get_samples_postprocess(cfg, "DAS")
    subdir = "SubPostFlows"
    mp.create_skeleton(subdir, None)
    mw.nf_make_post(
        "DAS", cfg, samples, conditions, subdir, "INFO", combinations=combos
    )
    produced = glob.glob(os.path.join(subdir, "*DAS_diego_subflow.nf"))
    assert produced, "generator produced no DIEGO subflow"
    text = Path(produced[0]).read_text()
    assert (
        "LOGS/${params.gSCOMBO}_${params.gCOMBO}_${params.gDASCOMPS}"
        "/DAS/diego/run_diego.log"
    ) in text
    assert "COMPSTR" not in text


def test_nf_piscem_skips_bam_postmapping_and_qc(workdir):
    cfg = _load("config_Test.json")
    cfg["WORKFLOWS"] = "QC,TRIMMING,MAPPING"
    cfg["MAPPING"]["TOOLS"] = {"piscem": "piscem"}
    cfg["MAPPING"]["Test"] = {"piscem": {"OPTIONS": {"INDEX": "", "MAP": ""}}}
    conditions = mp.get_conditions(cfg)
    _pre, sub, _post = mw.get_processes(cfg)
    combos = mw.get_combo(sub, cfg, conditions)
    samples = mp.get_samples(cfg)
    subdir = "SubFlows"
    mp.create_skeleton(subdir, None)
    mw.nf_make_sub(
        sub, cfg, samples, conditions, subdir, "INFO", combinations=combos
    )
    produced = glob.glob(os.path.join(subdir, "*_subflow.nf"))
    assert produced, "generator produced no Piscem subflow"
    workflow = Path(produced[0]).read_text().rsplit("workflow {", 1)[-1]
    assert "POSTMAPPING.out" not in workflow
    assert "QC_MAPPING(" not in workflow
    assert (
        "MULTIQC(QC_RAW.out.qc.concat(QC_TRIMMING.out.qc)"
        ".concat(MAPPING.out.logs).collect())"
    ) in workflow


# --------------------------------------------------------------------------- #
# Golden-file / snapshot tests: file-emitting generators
# --------------------------------------------------------------------------- #


def _generate(engine, workdir):
    cfg = _load("config_Test.json")
    conditions = mp.get_conditions(cfg)
    _pre, sub, _post = mw.get_processes(cfg)
    samples = mp.get_samples(cfg)
    combos = mw.get_combo(sub, cfg, conditions)
    if engine == "smk":
        subdir = "SubSnakes"
        mp.create_skeleton(subdir, None)
        mw.make_sub(sub, cfg, samples, conditions, subdir, "INFO", combinations=combos)
    else:
        subdir = "SubFlows"
        mp.create_skeleton(subdir, None)
        mw.nf_make_sub(
            sub, cfg, samples, conditions, subdir, "INFO", combinations=combos
        )
    return subdir


@pytest.mark.parametrize("engine", ["smk", "nf"])
def test_make_sub_golden(engine, workdir):
    subdir = _generate(engine, workdir)
    produced = sorted(
        os.path.relpath(p, workdir) for p in glob.glob(os.path.join(subdir, "*"))
    )
    assert produced, "generator produced no files"

    for rel in produced:
        raw = Path(workdir / rel).read_text()
        if rel.endswith(".json"):
            raw = json.dumps(json.loads(raw), indent=2, sort_keys=True)
        content = _normalize(raw, workdir)
        _check_golden(f"config_Test/{engine}/{os.path.basename(rel)}", content)


@pytest.mark.parametrize("engine", ["smk", "nf"])
def test_make_sub_file_set_stable(engine, workdir):
    subdir = _generate(engine, workdir)
    names = sorted(
        os.path.basename(p) for p in glob.glob(os.path.join(subdir, "*"))
    )
    _check_golden(f"config_Test/{engine}/_filelist.txt", "\n".join(names) + "\n")


def test_nf_rustqc_receives_deduplicated_bams(workdir):
    cfg = _load("config_Test.json")
    cfg["WORKFLOWS"] = "QC,TRIMMING,DEDUP,MAPPING"
    qc_options = cfg["QC"]["Test"]["fastqc"]
    cfg["QC"]["TOOLS"] = {"rustqc": "rustqc"}
    cfg["QC"]["Test"] = {"rustqc": qc_options}

    conditions = mp.get_conditions(cfg)
    _pre, subworkflows, _post = mw.get_processes(cfg)
    samples = mp.get_samples(cfg)
    combinations = mw.get_combo(subworkflows, cfg, conditions)
    subdir = "RustqcSubFlows"
    mp.create_skeleton(subdir, None)
    mw.nf_make_sub(
        subworkflows,
        cfg,
        samples,
        conditions,
        subdir,
        "INFO",
        combinations=combinations,
    )

    generated = "\n".join(
        path.read_text() for path in Path(subdir).glob("*.nf")
    )
    assert (
        "RUSTQC_MAPPING(POSTMAPPING.out.postmap.concat(POSTMAPPING.out.postmapuni)"
        ".concat(DEDUPBAM.out.dedup))"
    ) in generated


def test_nf_make_sub_passes_step_sequencing_overrides(workdir):
    cfg = _load("config_Test.json")
    cfg["QC"]["Test"]["fastqc"]["SEQUENCING"] = "single,rf"
    cfg["TRIMMING"]["SEQUENCING"] = "single,fr"
    cfg["TRIMMING"]["Test"]["bbduk"]["SEQUENCING"] = ""
    cfg["DEDUP"]["Test"]["picard"]["SEQUENCING"] = "single,rf"

    conditions = mp.get_conditions(cfg)
    _pre, subworkflows, _post = mw.get_processes(cfg)
    combinations = mw.get_combo(subworkflows, cfg, conditions)
    mp.create_skeleton("StepSequencingSubFlows", None)
    jobs = mw.nf_make_sub(
        subworkflows,
        cfg,
        mp.get_samples(cfg),
        conditions,
        "StepSequencingSubFlows",
        "INFO",
        combinations=combinations,
    )

    params = {
        os.path.basename(job[0]): job[3]
        for job in jobs
        if "fastqc-bbduk-picard" in os.path.basename(job[0])
    }
    assert params
    step_params = next(iter(params.values()))
    assert step_params["QCPAIRED"] == "single"
    assert step_params["QCSTRANDED"] == "rf"
    assert step_params["TRIMMINGPAIRED"] == "single"
    assert step_params["TRIMMINGSTRANDED"] == "fr"
    assert step_params["DEDUPPAIRED"] == "single"
    assert step_params["DEDUPSTRANDED"] == "rf"


# --------------------------------------------------------------------------- #
# Golden-file snapshot for the DTU postprocessing trio (drimseq/dexseq/spit).
# These share a large amount of scaffolding; the snapshot pins current output
# so the shared-base extraction can be reviewed as an explicit, minimal diff.
# --------------------------------------------------------------------------- #


def _generate_dtu(engine, workdir):
    cfg = _load("config_DTU_test.json")
    conditions = mp.get_conditions(cfg)
    _pre, sub, _post = mw.get_processes(cfg)
    combos = mw.get_combo(sub, cfg, conditions)
    samples = mp.get_samples_postprocess(cfg, "DTU")
    if engine == "smk":
        subdir = "SubPostSnakes"
        mp.create_skeleton(subdir, None)
        mw.make_post(
            "DTU", cfg, samples, conditions, subdir, "INFO", combinations=combos
        )
    else:
        subdir = "SubPostFlows"
        mp.create_skeleton(subdir, None)
        mw.nf_make_post(
            "DTU", cfg, samples, conditions, subdir, "INFO", combinations=combos
        )
    return subdir


def _sort_string_lists(node):
    """Sort lists whose members are all strings so set-derived orderings
    (e.g. postprocess SAMPLES) do not make the snapshot flaky."""
    if isinstance(node, dict):
        return {k: _sort_string_lists(v) for k, v in node.items()}
    if isinstance(node, list):
        items = [_sort_string_lists(v) for v in node]
        if all(isinstance(v, str) for v in items):
            return sorted(items)
        return items
    return node


@pytest.mark.parametrize("engine", ["smk", "nf"])
def test_make_post_dtu_golden(engine, workdir):
    subdir = _generate_dtu(engine, workdir)
    produced = sorted(
        os.path.relpath(p, workdir) for p in glob.glob(os.path.join(subdir, "*"))
    )
    assert produced, "generator produced no files"

    for rel in produced:
        raw = Path(workdir / rel).read_text()
        if rel.endswith(".json"):
            data = _sort_string_lists(json.loads(raw))
            raw = json.dumps(data, indent=2, sort_keys=True)
        content = _normalize(raw, workdir)
        _check_golden(f"config_DTU/{engine}/{os.path.basename(rel)}", content)


@pytest.mark.parametrize("engine", ["smk", "nf"])
def test_make_post_dtu_file_set_stable(engine, workdir):
    subdir = _generate_dtu(engine, workdir)
    names = sorted(
        os.path.basename(p) for p in glob.glob(os.path.join(subdir, "*"))
    )
    _check_golden(f"config_DTU/{engine}/_filelist.txt", "\n".join(names) + "\n")


def test_make_post_reuses_combinations_across_calls(workdir):
    """Regression: generating smk then nf postflows in one process must write
    files for both; get_combo used to return a one-shot iterator that the
    first call drained."""
    cfg = _load("config_DTU_test.json")
    conditions = mp.get_conditions(cfg)
    _pre, sub, _post = mw.get_processes(cfg)
    combos = mw.get_combo(sub, cfg, conditions)
    samples = mp.get_samples_postprocess(cfg, "DTU")

    mp.create_skeleton("SubPostSnakes", None)
    mw.make_post("DTU", cfg, samples, conditions, "SubPostSnakes", "INFO", combinations=combos)
    mp.create_skeleton("SubPostFlows", None)
    mw.nf_make_post("DTU", cfg, samples, conditions, "SubPostFlows", "INFO", combinations=combos)

    assert glob.glob(os.path.join("SubPostSnakes", "*")), "smk postflow empty"
    assert glob.glob(os.path.join("SubPostFlows", "*")), "nf postflow empty after smk call"


ILLUMINA_SAMPLES = [
    "381915_S1",
    "381915_S10",
    "381916_S2",
    "381917_S3",
    "381918_S4",
    "381919_S5",
]


@pytest.fixture
def illumina_workdir(tmp_path):
    """A temp working directory with real Illumina-style paired FASTQs.

    Only mate files carrying the ``_R1_001``/``_R2_001`` suffix exist, no
    canonical ``_R1.fastq.gz``/``_R2.fastq.gz`` names.
    """
    for item in os.listdir(DATA):
        if item == "FASTQ":
            continue
        os.symlink(DATA / item, tmp_path / item)
    fastq = tmp_path / "FASTQ" / "Test"
    fastq.mkdir(parents=True)
    for sample in ILLUMINA_SAMPLES:
        for read in ("R1", "R2"):
            (fastq / f"{sample}_{read}_001.fastq.gz").touch()
    old = os.getcwd()
    os.chdir(tmp_path)
    try:
        yield tmp_path
    finally:
        os.chdir(old)


@pytest.mark.parametrize("engine", ["smk", "nf"])
@pytest.mark.parametrize("subwork", ["DE", "DTU", "TRACKS"])
def test_postprocess_samples_strip_illumina_mate_suffix(
    subwork, engine, illumina_workdir
):
    """Regression: postprocess SAMPLES kept the _R1_001/_R2_001 mate suffix and
    reported every paired sample twice (12 IDs instead of 6)."""
    cfg = _load("config_Test.json")
    cfg["WORKFLOWS"] = "QC,TRIMMING,DEDUP,MAPPING,DE,DTU,TRACKS"
    cfg["SETTINGS"]["Test"]["SAMPLES"] = list(ILLUMINA_SAMPLES)
    cfg["SETTINGS"]["Test"]["GROUPS"] = ["T1", "T2", "T3", "T4", "T5", "T6"]
    expected = [f"Test/{sample}" for sample in ILLUMINA_SAMPLES]

    conditions = mp.get_conditions(cfg)
    _pre, sub, _post = mw.get_processes(cfg)
    combos = mw.get_combo(sub, cfg, conditions)
    samples = mp.get_samples_postprocess(cfg, subwork)
    assert sorted(samples) == sorted(expected)

    subdir = "SubPostSnakes" if engine == "smk" else "SubPostFlows"
    mp.create_skeleton(subdir, None)
    if engine == "smk":
        jobs = mw.make_post(
            subwork, cfg, samples, conditions, subdir, "INFO", combinations=combos
        )
    else:
        jobs = mw.nf_make_post(
            subwork, cfg, samples, conditions, subdir, "INFO", combinations=combos
        )

    confos = [Path(job[1]) for job in jobs]
    assert confos, "generator produced no subconfig"
    for confo in confos:
        assert sorted(json.loads(confo.read_text())["SAMPLES"]) == sorted(expected), confo

    for sample in ILLUMINA_SAMPLES:
        for read in ("R1", "R2"):
            mate = illumina_workdir / "FASTQ" / "Test" / f"{sample}_{read}_001.fastq.gz"
            assert mate.is_file(), mate


# --------------------------------------------------------------------------- #
# Golden-file snapshot for the FUSIONS postprocessing step (STAR-Fusion).
# --------------------------------------------------------------------------- #


def _generate_fusions(engine, workdir):
    cfg = _load("config_DTU_test.json")
    conditions = mp.get_conditions(cfg)
    _pre, sub, _post = mw.get_processes(cfg)
    combos = mw.get_combo(sub, cfg, conditions)
    samples = mp.get_samples_postprocess(cfg, "FUSIONS")
    if engine == "smk":
        subdir = "SubFusionSnakes"
        mp.create_skeleton(subdir, None)
        mw.make_post(
            "FUSIONS", cfg, samples, conditions, subdir, "INFO", combinations=combos
        )
    else:
        subdir = "SubFusionFlows"
        mp.create_skeleton(subdir, None)
        mw.nf_make_post(
            "FUSIONS", cfg, samples, conditions, subdir, "INFO", combinations=combos
        )
    return subdir


@pytest.mark.parametrize("engine", ["smk", "nf"])
def test_make_post_fusions_golden(engine, workdir):
    subdir = _generate_fusions(engine, workdir)
    produced = sorted(
        os.path.relpath(p, workdir) for p in glob.glob(os.path.join(subdir, "*"))
    )
    assert produced, "generator produced no files"

    for rel in produced:
        raw = Path(workdir / rel).read_text()
        if rel.endswith(".json"):
            data = _sort_string_lists(json.loads(raw))
            raw = json.dumps(data, indent=2, sort_keys=True)
        content = _normalize(raw, workdir)
        _check_golden(f"config_FUSIONS/{engine}/{os.path.basename(rel)}", content)


@pytest.mark.parametrize("engine", ["smk", "nf"])
def test_make_post_fusions_file_set_stable(engine, workdir):
    subdir = _generate_fusions(engine, workdir)
    names = sorted(
        os.path.basename(p) for p in glob.glob(os.path.join(subdir, "*"))
    )
    _check_golden(f"config_FUSIONS/{engine}/_filelist.txt", "\n".join(names) + "\n")



# --------------------------------------------------------------------------- #
# Cross-engine consistency: the Snakemake and Nextflow generators must expose
# the same tool matrix, env/bin selection and per-step OPTIONS for one config.
# The on-disk nesting differs between engines (smk nests OPTIONS under the tool
# name, nf flattens it), so comparisons are made on the semantic content.
# --------------------------------------------------------------------------- #

WORKFLOW_KEYS = ("QC", "TRIMMING", "MAPPING", "DEDUP", "COUNTING", "TRACKS")


def _combo_name(fname):
    for suffix in ("_subconfig.json", "_subsnake.smk", "_subflow.nf"):
        if fname.endswith(suffix):
            return fname[: -len(suffix)]
    return fname


def _subconfigs_by_combo(subdir):
    out = {}
    for path in glob.glob(os.path.join(subdir, "*_subconfig.json")):
        with open(path) as fh:
            out[_combo_name(os.path.basename(path))] = json.load(fh)
    return out


def _collect_options(node):
    """Recursively merge every OPTIONS dict found under a config subtree."""
    merged = {}
    if isinstance(node, dict):
        for key, val in node.items():
            if key == "OPTIONS" and isinstance(val, dict):
                merged.update(val)
            else:
                merged.update(_collect_options(val))
    return merged


@pytest.fixture
def both_engines(workdir):
    smk = _subconfigs_by_combo(_generate("smk", workdir))
    nf = _subconfigs_by_combo(_generate("nf", workdir))
    return smk, nf


def test_engines_expose_same_combos(both_engines):
    smk, nf = both_engines
    assert set(smk) == set(nf)
    assert smk, "no combos generated"


def test_engines_expose_same_env_and_bin(both_engines):
    smk, nf = both_engines
    for combo in smk:
        sc, nc = smk[combo], nf[combo]
        scalars = {
            k: v
            for k, v in sc.items()
            if k.endswith(("ENV", "BIN")) and isinstance(v, str)
        }
        for key, val in scalars.items():
            assert nc.get(key) == val, (
                f"{combo}: {key} drift smk={val!r} nf={nc.get(key)!r}"
            )


def test_engines_expose_same_options(both_engines):
    smk, nf = both_engines
    for combo in smk:
        sc, nc = smk[combo], nf[combo]
        for wf in WORKFLOW_KEYS:
            if wf not in sc and wf not in nc:
                continue
            smk_opts = _collect_options(sc.get(wf, {}))
            nf_opts = _collect_options(nc.get(wf, {}))
            assert smk_opts == nf_opts, (
                f"{combo}/{wf}: OPTIONS drift smk={smk_opts!r} nf={nf_opts!r}"
            )


def test_engines_expose_same_tools(both_engines):
    smk, nf = both_engines
    for combo in smk:
        sc, nc = smk[combo], nf[combo]
        for wf in WORKFLOW_KEYS:
            if wf not in sc and wf not in nc:
                continue
            assert sc.get(wf, {}).get("TOOLS") == nc.get(wf, {}).get("TOOLS"), (
                f"{combo}/{wf}: TOOLS drift"
            )


@pytest.mark.parametrize(
    "version,expected",
    [
        ("1.5.0", "1.5.0"),
        ("1.5.0+12.gabcdef", "1.5.0-12.gabcdef"),
        ("1.5.0+12.gabcdef.dirty", "1.5.0-12.gabcdef"),
        ("v2.0.0/a", "2.0.0-a"),
    ],
)
def test_normalize_container_version(version, expected):
    assert normalize_container_version(version) == expected


def test_normalize_container_version_rejects_empty_tag():
    with pytest.raises(ValueError):
        normalize_container_version("+++")


def test_oras_namespace_registry_and_version_rewrite(monkeypatch):
    line = 'container: "oras://jfallmann/monsda:"+MAPPERENV'
    monkeypatch.setattr(mw, "normalize_container_version", lambda version=None: "1.5.0")
    try:
        assert mw._rewrite_oras(line) == (
            'container: "oras://ghcr.io/jfallmann/monsda:"+MAPPERENV+"-1.5.0"'
        )
        mw.set_oras_namespace("myuser/myrepo")
        mw.set_oras_registry("registry.example.org")
        assert mw._rewrite_oras(line) == (
            'container: "oras://registry.example.org/myuser/myrepo:"+MAPPERENV+"-1.5.0"'
        )
    finally:
        mw.set_oras_namespace(mw.oras_default_namespace)
        mw.set_oras_registry("ghcr.io")


def test_oras_static_and_nextflow_containers_are_versioned(monkeypatch):
    monkeypatch.setattr(mw, "normalize_container_version", lambda version=None: "1.5.0")
    content = (
        'container: "oras://jfallmann/monsda:perl"\n'
        '    container "oras://jfallmann/monsda:"+"$QCENV"\n'
    )
    assert mw._rewrite_oras(content) == (
        'container: "oras://ghcr.io/jfallmann/monsda:perl"+"-1.5.0"\n'
        '    container "oras://ghcr.io/jfallmann/monsda:"+"$QCENV"+"-1.5.0"\n'
    )


def test_container_matrix_covers_every_environment():
    expected = sorted(
        path.stem for path in (REPO / "envs").glob("*.yaml")
        if path.stem not in SKIP_ENVIRONMENTS
    )
    assert environment_names(REPO) == expected


def test_container_build_disables_provenance_for_apptainer():
    workflow = yaml.safe_load(
        (REPO / ".github" / "workflows" / "build-containers.yml").read_text()
    )
    build_step = next(
        step
        for step in workflow["jobs"]["build"]["steps"]
        if step.get("uses") == "docker/build-push-action@v7"
    )
    assert build_step["with"]["provenance"] is False


def test_container_build_only_runs_for_releases_or_reusable_calls():
    workflow = yaml.safe_load(
        (REPO / ".github" / "workflows" / "build-containers.yml").read_text()
    )
    triggers = workflow.get("on", workflow.get(True))
    assert set(triggers) == {"release", "workflow_call"}
    assert triggers["release"]["types"] == ["published"]


def test_all_template_container_calls_are_versioned():
    calls = list()
    for pattern in ("*.smk", "*.nf"):
        for path in (REPO / "workflows").glob(pattern):
            calls.extend(
                line for line in path.read_text().splitlines()
                if "container" in line and "oras://" in line
            )
    assert calls
    rewritten = [mw._rewrite_oras(line) for line in calls]
    assert all(line.endswith('+"-VERSION"') for line in rewritten)


def test_static_container_environments_exist():
    static = set()
    pattern = re.compile(r"oras://jfallmann/monsda:([A-Za-z0-9_.-]+)\"")
    for suffix in ("*.smk", "*.nf"):
        for path in (REPO / "workflows").glob(suffix):
            static.update(pattern.findall(path.read_text()))
    available = set(environment_names(REPO))
    assert static
    assert static <= available


def test_rammap_is_provided_by_environment_without_runtime_install():
    environment = (REPO / "envs" / "rammap.yaml").read_text()
    assert "rammap =1.2.0" in environment
    for suffix in ("smk", "nf"):
        workflow = (REPO / "workflows" / f"rammap.{suffix}").read_text()
        assert "cargo install" not in workflow
        assert "RAMMAPBUILD" not in workflow
        assert "rammapbuild" not in workflow


def test_container_builds_make_command_scripts_executable():
    chmod_command = "find /opt/MONSDA/scripts -type f"
    definitions = [
        path
        for path in (REPO / "containers" / "apptainer").glob("*.def")
        if "scripts /opt/MONSDA" in path.read_text()
    ]
    assert definitions
    assert all(chmod_command in path.read_text() for path in definitions)
    assert chmod_command in (REPO / "docker" / "Dockerfile").read_text()


def test_scripts_called_without_interpreter_are_executable():
    called = re.compile(
        r"(?<![\w./])(?:\{params\.bins\}|\$\{?BINS\}?|\$\{params\.gBINS\})"
        r"/([\w/]+\.(?:py|pl|sh|R))"
    )
    interpreted = re.compile(r"(?:python3?|perl|Rscript|bash|sh)\s+$")
    direct = set()
    templates = sorted(
        set((REPO / "workflows").glob("*.smk")) | set((REPO / "workflows").glob("*.nf"))
    )
    for template in templates:
        text = template.read_text()
        for match in called.finditer(text):
            if not interpreted.search(text[: match.start()].rsplit("\n", 1)[-1]):
                direct.add(match.group(1))
    assert direct
    not_executable = sorted(
        name
        for name in direct
        if not os.access(REPO / "scripts" / name, os.X_OK)
    )
    assert not_executable == []


# --------------------------------------------------------------------------- #
# featureCounts paired-end counting: fragments, not reads.
# --------------------------------------------------------------------------- #

FC_PAIRED_FLAG = "-p --countReadPairs"

# workflow template -> (postprocessing step, generated subworkflow token, sites)
FC_TEMPLATES = {
    "countreads": ("COUNTING", "countreads", {"smk": 4, "nf": 1}),
    "deseq2_DE": ("DE", "deseq2", {"smk": 1, "nf": 1}),
    "edger_DE": ("DE", "edger", {"smk": 1, "nf": 1}),
    "dexseq_DEU": ("DEU", "dexseq", {"smk": 1, "nf": 1}),
    "edger_DEU": ("DEU", "edger", {"smk": 1, "nf": 1}),
    "edger_DAS": ("DAS", "edger", {"smk": 1, "nf": 1}),
    "diego_DAS": ("DAS", "diego", {"smk": 1, "nf": 1}),
}


def _fc_subworkflows(engine, layout, workdir, steps=("COUNTING", "DE", "DEU", "DAS")):
    """Generate the post subworkflows that carry a featureCounts branch."""
    cfg = _load("config_Test.json")
    cfg["SETTINGS"]["Test"]["SEQUENCING"] = layout
    cfg["DEU"]["TOOLS"]["dexseq"] = "Analysis/DEU/DEXSEQ.R"
    cfg["DEU"]["Test"]["dexseq"] = {
        "OPTIONS": {"QUANT": "-t exon -g gene_id -f", "DEU": ""}
    }
    cfg["DAS"]["TOOLS"]["diego"] = "Analysis/DAS/DIEGO.R"
    cfg["DAS"]["Test"]["diego"] = {
        "OPTIONS": {"QUANT": "-t exon -g gene_id -f", "DAS": ""}
    }
    conditions = mp.get_conditions(cfg)
    _pre, sub, _post = mw.get_processes(cfg)
    combos = mw.get_combo(sub, cfg, conditions)
    produced = {}
    for step in steps:
        subdir = f"Fc{engine}{layout}{step}"
        mp.create_skeleton(subdir, None)
        samples = mp.get_samples_postprocess(cfg, step)
        if engine == "smk":
            mw.make_post(
                step, cfg, samples, conditions, subdir, "INFO", combinations=combos
            )
            suffix = "subsnake.smk"
        else:
            mw.nf_make_post(
                step, cfg, samples, conditions, subdir, "INFO", combinations=combos
            )
            suffix = "subflow.nf"
        for conf in glob.glob(os.path.join(subdir, "*_subconfig.json")):
            sequencing = json.loads(Path(conf).read_text())["SETTINGS"]["Test"]
            assert sequencing["SEQUENCING"] == layout
        for path in glob.glob(os.path.join(subdir, f"*_{step}_*_{suffix}")):
            produced[os.path.basename(path)] = Path(path).read_text()
    return produced


def _fc_flags(text, engine, step, layout):
    """Resolve the featureCounts paired flag of every site in a subworkflow."""
    if engine == "smk":
        exprs = re.findall(
            r"paired\s*=\s*lambda x: (.+?) if %s_PAIRED == 'paired' else ''" % step,
            text,
        )
        return [eval(expr, {f"{step}_PAIRED": layout}) for expr in exprs]
    branches = re.findall(
        r"if \(params\.g%sPAIRED == 'paired'\)\{\s*pair = \"([^\"]*)\"\s*\}"
        r"\s*else\{\s*pair= \"([^\"]*)\"\s*\}" % step,
        text,
    )
    return [paired if layout == "paired" else single for paired, single in branches]


def _fc_rendered_command(text, engine):
    """The rendered featureCounts command line of the countreads subworkflow."""
    if engine == "smk":
        line = next(
            row.strip()
            for row in text.splitlines()
            if row.strip().startswith("shell:")
            and "{params.countb} -T {threads}" in row
        )
        return line.split("shell:", 1)[1].strip().strip('"')
    return next(
        row.strip()
        for row in text.splitlines()
        if "${params.gCOUNTBIN} -T ${task.cpus}" in row
    )


def _fc_run(engine, text, layout, fc, workdir):
    """Run the rendered countreads command with the real featureCounts."""
    run = workdir / f"fc_run_{engine}_{layout}"
    run.mkdir()
    (run / "LOGS").mkdir()
    anno = str(workdir / "fc_anno.gtf.gz")
    reads = str(workdir / ("fc_frag.sam" if layout == "paired" else "fc_single.sam"))
    flag = FC_PAIRED_FLAG if layout == "paired" else ""
    command = _fc_rendered_command(text, engine)
    env = None
    if engine == "smk":
        for key, value in {
            "{params.countb}": fc,
            "{threads}": "1",
            "{params.cpara}": "-t exon -g gene_id",
            "{params.paired}": flag,
            "{params.stranded}": "",
            "{params.anno}": anno,
            "{output.t}": "tmpcts",
            "{input.s}": reads,
            "{log}": "LOGS/count.log",
            "{output.cts}": "counts.gz",
            "{params.sortmem}": "1",
        }.items():
            command = command.replace(key, value)
    else:
        command = (
            command.replace("${params.gCOUNTBIN}", fc)
            .replace("${params.gCOUNTPARAMS}", "-t exon -g gene_id")
            .replace("${params.gCOUNTMAP}", "")
            .replace("${task.cpus}", "1")
        )
        env = dict(
            os.environ,
            pair=flag,
            stranded="",
            anno=anno,
            reads=reads,
            ol="count.log",
            oc="counts.gz",
            os="counts.summary",
            sortmem="1",
        )
    result = subprocess.run(
        ["bash", "-c", command], cwd=run, env=env, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr
    rows = gzip.open(run / "counts.gz", "rt").read().splitlines()[2:]
    return {row.split("\t")[0]: int(row.split("\t")[-1]) for row in rows}


@pytest.mark.parametrize("engine", ["smk", "nf"])
@pytest.mark.parametrize("layout", ["paired", "single"])
def test_featurecounts_paired_flag_counts_fragments(engine, layout, workdir):
    produced = _fc_subworkflows(engine, layout, workdir)
    expected = FC_PAIRED_FLAG if layout == "paired" else ""
    for template, (step, token, sites) in FC_TEMPLATES.items():
        name = next(n for n in produced if f"_{step}_{token}_sub" in n)
        flags = _fc_flags(produced[name], engine, step, layout)
        assert len(flags) == sites[engine], (template, flags)
        assert set(flags) == {expected}, (template, flags)


@pytest.mark.parametrize("engine", ["smk", "nf"])
def test_featurecounts_rendered_command_counts_fragments(engine, workdir):
    fc = shutil.which("featureCounts")
    if not fc:
        pytest.skip("featureCounts is not on PATH")

    anno = workdir / "fc_anno.gtf"
    anno.write_text(
        'chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    )
    subprocess.run(["gzip", "-kf", str(anno)], check=True)
    seq, qual = "A" * 50, "I" * 50
    header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n"
    (workdir / "fc_frag.sam").write_text(
        header
        + f"r1\t99\tchr1\t120\t60\t50M\t=\t160\t90\t{seq}\t{qual}\n"
        + f"r1\t147\tchr1\t160\t60\t50M\t=\t120\t-90\t{seq}\t{qual}\n"
    )
    (workdir / "fc_single.sam").write_text(
        header
        + f"s1\t0\tchr1\t120\t60\t50M\t*\t0\t0\t{seq}\t{qual}\n"
        + f"s2\t0\tchr1\t160\t60\t50M\t*\t0\t0\t{seq}\t{qual}\n"
    )

    expected = {"paired": 1, "single": 2}
    for layout in ("paired", "single"):
        produced = _fc_subworkflows(engine, layout, workdir, steps=("COUNTING",))
        name = next(n for n in produced if "_COUNTING_countreads_sub" in n)
        counts = _fc_run(engine, produced[name], layout, fc, workdir)
        assert counts["g1"] == expected[layout], (engine, layout, counts)
