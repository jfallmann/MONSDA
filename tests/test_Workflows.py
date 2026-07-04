import glob
import json
import os
from pathlib import Path

import pytest

pytest.importorskip("snakemake")

from snakemake.common.configfile import load_configfile  # noqa: E402

import MONSDA.Params as mp  # noqa: E402
import MONSDA.Workflows as mw  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "tests" / "data"
GOLDEN = REPO / "tests" / "golden"

UPDATE = bool(os.environ.get("UPDATE_GOLDEN"))


@pytest.fixture(autouse=True)
def _repo_template_paths():
    """Point the generators at the in-repo workflow/env templates.

    In a dev checkout MONSDA is not pip-installed, so the module-level
    ``workflowpath``/``envpath`` (which assume a site-packages layout) do not
    resolve. Redirect them to the repository copies for the duration of a test.
    """
    old_wf, old_env = mw.workflowpath, mw.envpath
    mw.workflowpath = str(REPO / "workflows")
    mw.envpath = str(REPO / "envs") + os.sep
    yield
    mw.workflowpath, mw.envpath = old_wf, old_env


def _load(name):
    return load_configfile(str(DATA / name))


def _normalize(text, workdir):
    return text.replace(str(workdir), "<WORKDIR>").replace(str(REPO), "<REPO>")


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


def test_oras_namespace_and_registry_rewrite():
    line = 'container: "oras://jfallmann/monsda:"+MAPPERENV'
    try:
        assert mw._rewrite_oras(line) == (
            'container: "oras://docker.io/jfallmann/monsda:"+MAPPERENV'
        )
        mw.set_oras_namespace("myuser/myrepo")
        mw.set_oras_registry("ghcr.io")
        assert mw._rewrite_oras(line) == (
            'container: "oras://ghcr.io/myuser/myrepo:"+MAPPERENV'
        )
    finally:
        mw.set_oras_namespace(mw.oras_default_namespace)
        mw.set_oras_registry("docker.io")
