import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("snakemake")

from snakemake.common.configfile import load_configfile  # noqa: E402

import MONSDA.Params as mp  # noqa: E402
import MONSDA.Workflows as mw  # noqa: E402
from MONSDA.PostDE import load_postde_config, prepare_postde  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "tests" / "data"


@pytest.fixture(autouse=True)
def _repo_template_paths(monkeypatch):
    old_wf, old_env, old_bin = mw.workflowpath, mw.envpath, mw.binpath
    mw.workflowpath = str(REPO / "workflows")
    mw.envpath = str(REPO / "envs") + os.sep
    mw.binpath = str(REPO / "scripts")
    monkeypatch.setattr(
        mw, "normalize_container_version", lambda version=None: "VERSION"
    )
    yield
    mw.workflowpath, mw.envpath, mw.binpath = old_wf, old_env, old_bin


def _load(name):
    return load_configfile(str(DATA / name))


@pytest.fixture
def workdir(tmp_path):
    for item in os.listdir(DATA):
        os.symlink(DATA / item, tmp_path / item)
    old = os.getcwd()
    os.chdir(tmp_path)
    try:
        yield tmp_path
    finally:
        os.chdir(old)


def _analysis_config(tmp_path, analysis, name="postde_analysis.json"):
    path = tmp_path / name
    path.write_text(json.dumps(analysis))
    return path


def _postde_cfg(tmp_path, analysis, workflows="DE", tools=None):
    cfg = _load("config_Test.json")
    cfg["WORKFLOWS"] = workflows
    if tools is not None:
        cfg["DE"]["TOOLS"] = tools
    cfg["POSTDE"] = {
        "enabled": True,
        "config": str(_analysis_config(tmp_path, analysis)),
    }
    return cfg


def _generate_de(engine, workdir, cfg):
    conditions = mp.get_conditions(cfg)
    samples = mp.get_samples_postprocess(cfg, "DE")
    if engine == "smk":
        subdir = "SubSnakes"
        mp.create_skeleton(subdir, None)
        mw.make_post("DE", cfg, samples, conditions, subdir, "INFO")
    else:
        subdir = "SubFlows"
        mp.create_skeleton(subdir, None)
        mw.nf_make_post("DE", cfg, samples, conditions, subdir, "INFO")
    return subdir


def _extract_nf_script(subflow, process_name):
    match = re.search(r"process\s+" + process_name + r"\s*\{", subflow)
    assert match, f"process {process_name} not found in generated subflow"
    start = match.end()
    depth = 1
    i = start
    while i < len(subflow) and depth > 0:
        if subflow[i] == "{":
            depth += 1
        elif subflow[i] == "}":
            depth -= 1
        i += 1
    body = subflow[start:i]
    block = re.search(r'script:\s*"""\n(.*?)\n\s*"""', body, re.S)
    assert block, f"script block not found in process {process_name}"
    return block.group(1)


def _write_stub_bundle(rscript, path):
    rcode = (
        'bundle <- list(schema_version = 1L, engine = "deseq2", '
        'contrasts = list(list(id = "T1-VS-T2")))\n'
        "saveRDS(bundle, commandArgs(trailingOnly = TRUE)[1])\n"
    )
    rfile = path.with_suffix(".R")
    rfile.write_text(rcode)
    subprocess.run(
        [rscript, str(rfile), str(path)], check=True, capture_output=True
    )


class TestLoadPostdeConfig:
    def test_none_returns_none(self):
        assert load_postde_config(None) is None

    def test_disabled_returns_section(self):
        section = {"enabled": False, "config": ""}
        assert load_postde_config(section) == section

    def test_missing_enabled_raises(self):
        with pytest.raises(ValueError, match="enabled"):
            load_postde_config({"config": ""})

    def test_enabled_nonbool_raises(self):
        with pytest.raises(ValueError, match="boolean"):
            load_postde_config({"enabled": "yes", "config": ""})

    def test_unknown_section_field_raises(self):
        with pytest.raises(ValueError, match="unknown fields"):
            load_postde_config({"enabled": False, "config": "", "enabeld": True})

    def test_enabled_missing_config_raises(self):
        with pytest.raises(ValueError, match="config"):
            load_postde_config({"enabled": True})

    def test_config_not_found_raises(self):
        with pytest.raises(ValueError, match="not found"):
            load_postde_config({"enabled": True, "config": "nope.json"})

    def test_config_invalid_json_raises(self, tmp_path):
        path = tmp_path / "bad.json"
        path.write_text("not json")
        with pytest.raises(ValueError, match="JSON"):
            load_postde_config({"enabled": True, "config": str(path)})

    def test_unknown_analysis_field_raises(self, tmp_path):
        path = _analysis_config(tmp_path, {"gprofilerr": {"enabled": True}})
        with pytest.raises(ValueError, match="unknown top-level"):
            load_postde_config({"enabled": True, "config": str(path)})

    def test_plots_enabled_raises(self, tmp_path):
        path = _analysis_config(tmp_path, {"plots": {"enabled": True}})
        with pytest.raises(ValueError, match="not available"):
            load_postde_config({"enabled": True, "config": str(path)})

    def test_gprofiler_requires_organism(self, tmp_path):
        path = _analysis_config(
            tmp_path, {"gprofiler": {"enabled": True, "allow_network": True}}
        )
        with pytest.raises(ValueError, match="organism"):
            load_postde_config({"enabled": True, "config": str(path)})

    def test_gprofiler_requires_allow_network(self, tmp_path):
        path = _analysis_config(
            tmp_path, {"gprofiler": {"enabled": True, "organism": "hsapiens"}}
        )
        with pytest.raises(ValueError, match="allow_network"):
            load_postde_config({"enabled": True, "config": str(path)})

    def test_decoupler_requires_network_or_allow_network(self, tmp_path):
        path = _analysis_config(tmp_path, {"decoupler": {"enabled": True}})
        with pytest.raises(ValueError, match="allow_network"):
            load_postde_config({"enabled": True, "config": str(path)})

    def test_missing_resource_raises(self, tmp_path):
        path = _analysis_config(
            tmp_path,
            {"clusterprofiler": {"enabled": True, "term2gene": "missing.tsv"}},
        )
        with pytest.raises(ValueError, match="resource file not found"):
            load_postde_config({"enabled": True, "config": str(path)})

    def test_valid_returns_abs_config(self, tmp_path):
        (tmp_path / "term2gene.tsv").write_text("term\tgene\nGO:1\tG1\n")
        analysis = {
            "clusterprofiler": {"enabled": True, "term2gene": "term2gene.tsv"}
        }
        path = _analysis_config(tmp_path, analysis)
        out = load_postde_config({"enabled": True, "config": str(path)})
        assert out["enabled"] is True
        assert out["config"] == str(path.resolve())

    def test_no_de_in_workflows_raises(self, tmp_path):
        cfg = _postde_cfg(tmp_path, {}, workflows="QC")
        with pytest.raises(ValueError, match="DE in WORKFLOWS"):
            load_postde_config(cfg["POSTDE"], config=cfg)

    def test_unsupported_de_tool_raises(self, tmp_path):
        cfg = _postde_cfg(
            tmp_path, {}, tools={"diego": "Analysis/DE/DIEGO.R"}
        )
        with pytest.raises(ValueError, match="does not support"):
            load_postde_config(cfg["POSTDE"], config=cfg)


class TestPreparePostde:
    def test_disabled_returns_none(self, tmp_path):
        cfg = _load("config_Test.json")
        cfg["POSTDE"] = {"enabled": False, "config": ""}
        assert prepare_postde(cfg, "SubSnakes") is None
        assert not (tmp_path / "SubSnakes").exists()

    def test_no_postde_returns_none(self, tmp_path):
        cfg = _load("config_Test.json")
        assert prepare_postde(cfg, "SubSnakes") is None

    def test_stages_config_and_resources(self, tmp_path):
        (tmp_path / "term2gene.tsv").write_text("term\tgene\nGO:1\tG1\n")
        (tmp_path / "term2name.tsv").write_text("term\tname\nGO:1\tName1\n")
        analysis = {
            "padj": 0.05,
            "clusterprofiler": {
                "enabled": True,
                "term2gene": "term2gene.tsv",
                "term2name": "term2name.tsv",
            },
        }
        cfg = _postde_cfg(tmp_path, analysis)
        staged = prepare_postde(cfg, "SubSnakes")
        assert staged is not None
        assert staged.startswith("SubSnakes" + os.sep + "POSTDE_")
        assert os.path.isfile(os.path.join(staged, "config.json"))
        assert os.path.isfile(os.path.join(staged, "term2gene.tsv"))
        assert os.path.isfile(os.path.join(staged, "term2name.tsv"))
        staged_cfg = json.load(open(os.path.join(staged, "config.json")))
        assert staged_cfg["clusterprofiler"]["term2gene"] == "term2gene.tsv"
        assert staged_cfg["clusterprofiler"]["term2name"] == "term2name.tsv"

    def test_relative_resource_paths(self, tmp_path):
        resdir = tmp_path / "res"
        resdir.mkdir()
        (resdir / "background.tsv").write_text("G1\n")
        analysis = {
            "gprofiler": {
                "enabled": True,
                "organism": "hsapiens",
                "allow_network": True,
                "background": "res/background.tsv",
            }
        }
        cfg = _postde_cfg(tmp_path, analysis)
        staged = prepare_postde(cfg, "SubSnakes")
        assert os.path.isfile(os.path.join(staged, "background.tsv"))

    def test_duplicate_basenames(self, tmp_path):
        (tmp_path / "background.tsv").write_text("G1\n")
        resdir = tmp_path / "res"
        resdir.mkdir()
        (resdir / "background.tsv").write_text("G2\n")
        analysis = {
            "gprofiler": {
                "enabled": True,
                "organism": "hsapiens",
                "allow_network": True,
                "background": "background.tsv",
            },
            "clusterprofiler": {
                "enabled": True,
                "term2gene": "res/background.tsv",
            },
        }
        cfg = _postde_cfg(tmp_path, analysis)
        staged = prepare_postde(cfg, "SubSnakes")
        names = sorted(os.listdir(staged))
        assert "background.tsv" in names
        assert "term2gene_background.tsv" in names

    def test_content_change_new_dir(self, tmp_path):
        (tmp_path / "term2gene.tsv").write_text("term\tgene\nGO:1\tG1\n")
        analysis = {
            "clusterprofiler": {"enabled": True, "term2gene": "term2gene.tsv"}
        }
        cfg = _postde_cfg(tmp_path, analysis)
        staged1 = prepare_postde(cfg, "SubSnakes")
        (tmp_path / "term2gene.tsv").write_text("term\tgene\nGO:1\tG2\n")
        staged2 = prepare_postde(cfg, "SubSnakes")
        assert staged1 != staged2

    def test_reuse_existing_dir(self, tmp_path):
        (tmp_path / "term2gene.tsv").write_text("term\tgene\nGO:1\tG1\n")
        analysis = {
            "clusterprofiler": {"enabled": True, "term2gene": "term2gene.tsv"}
        }
        cfg = _postde_cfg(tmp_path, analysis)
        staged1 = prepare_postde(cfg, "SubSnakes")
        mtime1 = os.path.getmtime(os.path.join(staged1, "config.json"))
        staged2 = prepare_postde(cfg, "SubSnakes")
        assert staged1 == staged2
        assert os.path.getmtime(os.path.join(staged2, "config.json")) == mtime1


class TestWorkflowIntegration:
    def test_make_post_smk(self, workdir, tmp_path):
        (tmp_path / "term2gene.tsv").write_text("term\tgene\nGO:1\tG1\n")
        analysis = {
            "clusterprofiler": {"enabled": True, "term2gene": "term2gene.tsv"}
        }
        cfg = _postde_cfg(tmp_path, analysis)
        subdir = _generate_de("smk", workdir, cfg)
        confo = json.load(
            open(workdir / subdir / "allconditions_DE_deseq2_subconfig.json")
        )
        assert confo["POSTDE"]["enabled"] is True
        assert confo["POSTDE"]["inputs"].startswith("SubSnakes" + os.sep + "POSTDE_")
        smk = (workdir / subdir / "allconditions_DE_deseq2_subsnake.smk").read_text()
        assert "rule postde" in smk
        assert "MONSDA_POSTDE={postde_flag}" in smk
        assert "DE_deseq2_{scombo}_postde.rds" in smk
        assert "POSTDE/{combo}/manifest.json" in smk
        smk_e = (workdir / subdir / "allconditions_DE_edger_subsnake.smk").read_text()
        assert "rule postde" in smk_e
        assert "DE_edger_{scombo}_postde.rds" in smk_e

    def test_make_post_nf(self, workdir, tmp_path):
        (tmp_path / "term2gene.tsv").write_text("term\tgene\nGO:1\tG1\n")
        analysis = {
            "clusterprofiler": {"enabled": True, "term2gene": "term2gene.tsv"}
        }
        cfg = _postde_cfg(tmp_path, analysis)
        subdir = _generate_de("nf", workdir, cfg)
        confo = json.load(
            open(workdir / subdir / "allconditions_DE_deseq2_subconfig.json")
        )
        assert confo["POSTDE"]["enabled"] is True
        assert confo["POSTDE"]["inputs"].startswith("SubFlows" + os.sep + "POSTDE_")
        nf = (workdir / subdir / "allconditions_DE_deseq2_subflow.nf").read_text()
        assert "process postde" in nf
        assert "MONSDA_POSTDE=${params.gPOSTDE_FLAG}" in nf
        assert "emit: bundle, optional: true" in nf
        assert "postde(run_deseq2.out.bundle" in nf
        nf_e = (workdir / subdir / "allconditions_DE_edger_subflow.nf").read_text()
        assert "process postde" in nf_e
        assert "postde(run_edger.out.bundle" in nf_e

    def test_disabled_no_postde(self, workdir, tmp_path):
        cfg = _load("config_Test.json")
        cfg["POSTDE"] = {"enabled": False, "config": ""}
        subdir = _generate_de("smk", workdir, cfg)
        confo = json.load(
            open(workdir / subdir / "allconditions_DE_deseq2_subconfig.json")
        )
        assert "POSTDE" not in confo
        smk = (workdir / subdir / "allconditions_DE_deseq2_subsnake.smk").read_text()
        assert "rule postde" not in smk
        assert "MONSDA_POSTDE=0" in smk

    def test_snakemake_dag_dryrun(self, workdir, tmp_path):
        (tmp_path / "term2gene.tsv").write_text("term\tgene\nGO:1\tG1\n")
        analysis = {
            "clusterprofiler": {"enabled": True, "term2gene": "term2gene.tsv"}
        }
        cfg = _postde_cfg(tmp_path, analysis)
        subdir = _generate_de("smk", workdir, cfg)
        confo = json.load(
            open(workdir / subdir / "allconditions_DE_deseq2_subconfig.json")
        )
        bundle = workdir / "DE" / "deseq2_DE" / "DE_deseq2__postde.rds"
        bundle.parent.mkdir(parents=True, exist_ok=True)
        bundle.write_bytes(b"fake")
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "snakemake",
                "-n",
                "-s",
                str(workdir / subdir / "allconditions_DE_deseq2_subsnake.smk"),
                "--configfile",
                str(workdir / subdir / "allconditions_DE_deseq2_subconfig.json"),
                "--directory",
                str(workdir),
                "POSTDE/deseq2_DE/manifest.json",
            ],
            capture_output=True,
            text=True,
            env={**os.environ, "PYTHONPATH": str(REPO)},
        )
        assert result.returncode == 0, result.stderr[-2000:]
        assert "postde" in result.stdout


class TestNfProcessBody:
    def test_postde_process_body_runs(self, workdir, tmp_path):
        rscript = shutil.which("Rscript")
        if not rscript:
            pytest.skip("Rscript not found on PATH")
        analysis = {"padj": 0.05, "lfc": 1}
        cfg = _postde_cfg(tmp_path, analysis)
        subdir = _generate_de("nf", workdir, cfg)
        confo = json.load(
            open(workdir / subdir / "allconditions_DE_deseq2_subconfig.json")
        )
        staged = confo["POSTDE"]["inputs"]
        subflow = (workdir / subdir / "allconditions_DE_deseq2_subflow.nf").read_text()
        script = _extract_nf_script(subflow, "postde")
        bundle = workdir / "DE" / "deseq2_DE" / "DE_deseq2__postde.rds"
        bundle.parent.mkdir(parents=True, exist_ok=True)
        _write_stub_bundle(rscript, bundle)
        outdir = workdir / "postde_out"
        script = script.replace("${params.gBINS}", str(REPO / "scripts"))
        script = script.replace("postde_inputs", str(workdir / staged))
        script = script.replace("../bundle", str(bundle))
        script = script.replace("../postde", str(outdir))
        result = subprocess.run(
            script, shell=True, cwd=str(workdir), capture_output=True, text=True
        )
        assert result.returncode == 0, result.stderr[-2000:]
        assert (outdir / "manifest.json").exists()


class TestCliSave:
    def test_cli_save_generates_postde(self, tmp_path):
        for item in ("FASTQ", "GENOME"):
            os.symlink(DATA / item, tmp_path / item)
        (tmp_path / "sitecustomize.py").write_text(
            "import MONSDA.Workflows as mw\n"
            f"mw.workflowpath = {str(REPO / 'workflows')!r}\n"
            f"mw.envpath = {str(REPO / 'envs') + os.sep!r}\n"
            f"mw.binpath = {str(REPO / 'scripts')!r}\n"
        )
        from MONSDA import _version

        version = _version.get_versions()["version"]
        cfg = _load("config_Test.json")
        cfg["VERSION"] = version
        cfg["WORKFLOWS"] = "DE"
        cfg["POSTDE"] = {
            "enabled": True,
            "config": str(_analysis_config(tmp_path, {"padj": 0.05, "lfc": 1})),
        }
        (tmp_path / "config_Test.json").write_text(json.dumps(cfg))
        env = {**os.environ, "PYTHONPATH": str(REPO) + os.pathsep + str(tmp_path)}
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "MONSDA.RunMONSDA",
                "--save",
                "-c",
                "config_Test.json",
            ],
            cwd=str(tmp_path),
            env=env,
            capture_output=True,
            text=True,
            timeout=300,
        )
        assert result.returncode == 0, result.stderr[-2000:]
        subdir = tmp_path / "SubSnakes"
        confo = json.load(open(subdir / "allconditions_DE_deseq2_subconfig.json"))
        assert confo["POSTDE"]["enabled"] is True
        staged = confo["POSTDE"]["inputs"]
        assert (tmp_path / staged / "config.json").exists()
        commands = (tmp_path / "JOBS" / "MONSDA.commands").read_text()
        assert "allconditions_DE_deseq2_subsnake.smk" in commands
