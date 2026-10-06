import json
import os
import re
import shlex
import shutil
import string
import subprocess
import sys
import types
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


@pytest.fixture(autouse=True)
def _sandbox_chdir(tmp_path, monkeypatch):
    """All staging tests run inside tmp_path so the repo tree is never touched."""
    monkeypatch.chdir(tmp_path)


def _load(name):
    return load_configfile(str(DATA / name))


@pytest.fixture
def workdir(tmp_path):
    for item in os.listdir(DATA):
        src = DATA / item
        if src.is_dir():
            os.symlink(src, tmp_path / item)
        else:
            shutil.copy2(src, tmp_path / item)
    return tmp_path


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


def _write_real_bundle(rscript, outdir, engine):
    """Create a tiny valid postde bundle via export.R: 50 genes x 4 samples
    (2 conditions, 2 reps), one contrast A-VS-B. Returns the bundle path."""
    export_r = REPO / "scripts" / "Analysis" / "PostDE" / "export.R"
    rcode = (
        "args <- commandArgs(trailingOnly = TRUE)\n"
        "export_r <- args[1]\n"
        "outdir <- args[2]\n"
        "engine <- args[3]\n"
        "postde_lib <- export_r\n"
        "source(export_r)\n"
        "set.seed(42)\n"
        "genes <- paste0('G', seq_len(50))\n"
        "samples <- c('S1', 'S2', 'S3', 'S4')\n"
        "condition <- factor(c('A', 'A', 'B', 'B'), levels = c('A', 'B'))\n"
        "metadata <- data.frame(row.names = samples, condition = condition)\n"
        "counts <- matrix(rnbinom(50 * 4, mu = 100, size = 10), nrow = 50, "
        "dimnames = list(genes, samples))\n"
        "expression <- matrix(rnorm(50 * 4), nrow = 50, "
        "dimnames = list(genes, samples))\n"
        "results <- data.frame(gene_id = genes, logFC = rnorm(50), "
        "pvalue = runif(50), padj = runif(50), stat = rnorm(50), "
        "mean = runif(50, 1, 100), row.names = genes, "
        "stringsAsFactors = FALSE)\n"
        "formula <- ~condition\n"
        "postde_capture(engine = engine, id = 'A-VS-B', A = 'A', B = 'B', "
        "normalized = FALSE, metadata = metadata, counts = counts, "
        "expression = expression, formula = formula, results = results, "
        "mean_scale = 'baseMean')\n"
        "postde_write(outdir, combi = '', engine = engine)\n"
    )
    rfile = outdir / ("make_bundle_" + engine + ".R")
    rfile.write_text(rcode)
    subprocess.run(
        [rscript, str(rfile), str(export_r), str(outdir), engine],
        check=True,
        capture_output=True,
    )
    return outdir / ("DE_" + engine + "__postde.rds")


def _extract_smk_rule_body(smk, rule_name):
    match = re.search(r"^\s*rule\s+" + rule_name + r"\s*:\s*$", smk, re.M)
    assert match, f"rule {rule_name} not found in generated snakefile"
    rest = smk[match.end():]
    body = []
    for ln in rest.split("\n"):
        if body and re.match(
            r"^\s*(rule|if|include|configfile|import|from|onsuccess|onerror)\b", ln
        ):
            break
        body.append(ln)
    return "\n".join(body)


def _smk_params_shell(body):
    params_m = re.search(r"params:\s*(.*?)(?=\n\s*shell:)", body, re.S)
    assert params_m, "params block not found in rule body"
    shell_m = re.search(r'shell:\s*"(.*?)"', body, re.S)
    assert shell_m, "shell block not found in rule body"
    lines = params_m.group(1).split("\n")
    code = "\n".join(
        [lines[0].strip().rstrip(",")]
        + [ln.strip().rstrip(",") for ln in lines[1:] if ln.strip()]
    )
    return code, shell_m.group(1)


def _format_shell(template, **fields):
    class _QFmt(string.Formatter):
        def format_field(self, value, spec):
            if spec == "q":
                return shlex.quote(str(value))
            return super().format_field(value, spec)

    return _QFmt().format(template, **fields)


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
        assert staged.startswith(os.path.abspath("SubSnakes") + os.sep + "POSTDE_")
        assert os.path.isfile(os.path.join(staged, "config.json"))
        assert os.path.isfile(os.path.join(staged, "clusterprofiler_term2gene.tsv"))
        assert os.path.isfile(os.path.join(staged, "clusterprofiler_term2name.tsv"))
        staged_cfg = json.load(open(os.path.join(staged, "config.json")))
        assert (
            staged_cfg["clusterprofiler"]["term2gene"]
            == "clusterprofiler_term2gene.tsv"
        )
        assert (
            staged_cfg["clusterprofiler"]["term2name"]
            == "clusterprofiler_term2name.tsv"
        )

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
        assert os.path.isfile(os.path.join(staged, "gprofiler_background.tsv"))

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
        assert "gprofiler_background.tsv" in names
        assert "clusterprofiler_term2gene.tsv" in names
        assert "background.tsv" not in names

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

    def test_cache_config_mutation_repaired(self, tmp_path):
        (tmp_path / "term2gene.tsv").write_text("term\tgene\nGO:1\tG1\n")
        analysis = {
            "padj": 0.05,
            "clusterprofiler": {"enabled": True, "term2gene": "term2gene.tsv"},
        }
        cfg = _postde_cfg(tmp_path, analysis)
        staged1 = prepare_postde(cfg, "SubSnakes")
        cfg_path = os.path.join(staged1, "config.json")
        original = open(cfg_path).read()
        # same-length in-place mutation of the staged config content
        mutated = original.replace('"padj": 0.05', '"padj": 0.06')
        assert len(mutated) == len(original)
        with open(cfg_path, "w") as fh:
            fh.write(mutated)
        # preserve mtime to defeat stat-based caching
        st = os.stat(cfg_path)
        os.utime(cfg_path, (st.st_atime, st.st_mtime))
        staged2 = prepare_postde(cfg, "SubSnakes")
        assert staged1 == staged2
        assert open(cfg_path).read() == original


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
        assert confo["POSTDE"]["inputs"].startswith(
            os.path.abspath("SubSnakes") + os.sep + "POSTDE_"
        )
        smk = (workdir / subdir / "allconditions_DE_deseq2_subsnake.smk").read_text()
        assert "rule postde" in smk
        assert "MONSDA_POSTDE={postde_flag}" in smk
        assert "DE_deseq2_{scombo}_postde.rds" in smk
        assert "POSTDE/{combo}" in smk
        assert "manifest.json" in smk
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
        assert confo["POSTDE"]["inputs"].startswith(
            os.path.abspath("SubFlows") + os.sep + "POSTDE_"
        )
        nf = (workdir / subdir / "allconditions_DE_deseq2_subflow.nf").read_text()
        assert "process postde" in nf
        assert "MONSDA_POSTDE=${params.gPOSTDE_FLAG}" in nf
        assert "emit: bundle, optional: !params.gPOSTDE_ENABLED" in nf
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
        assert "MONSDA_POSTDE={postde_flag}" in smk
        base = [
            sys.executable,
            "-m",
            "snakemake",
            "-s",
            str(workdir / subdir / "allconditions_DE_deseq2_subsnake.smk"),
            "--configfile",
            str(workdir / subdir / "allconditions_DE_deseq2_subconfig.json"),
            "--directory",
            str(workdir),
        ]
        env = {**os.environ, "PYTHONPATH": str(REPO)}
        listed = subprocess.run(
            base + ["--list-rules"], capture_output=True, text=True, env=env
        )
        assert listed.returncode == 0, listed.stderr[-2000:]
        assert "postde" not in listed.stdout
        dag = subprocess.run(
            base + ["-n", "POSTDE/deseq2_DE"], capture_output=True, text=True, env=env
        )
        assert dag.returncode != 0
        assert "MissingRuleException" in dag.stderr or "No rule to produce" in dag.stderr

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
                "POSTDE/deseq2_DE",
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
        # simulate nextflow rendering of the escaped shell variables
        script = script.replace("\\$", "$")
        script = script.replace("${params.gBINS}", str(REPO / "scripts"))
        _write_stub_bundle(rscript, workdir / "bundle")
        os.symlink(workdir / staged, workdir / "postde_inputs")
        result = subprocess.run(
            script, shell=True, cwd=str(workdir), capture_output=True, text=True
        )
        assert result.returncode == 0, result.stderr[-2000:]
        assert (workdir / "postde" / "manifest.json").exists()
        assert (workdir / "log").exists()


class TestPostdeDefaults:
    def test_gsva_and_global_minmax_defaults(self, tmp_path):
        """Effective GSVA/global min/max defaults must match the R
        common/enrichment defaults (10/500), not arbitrary values."""
        rscript = shutil.which("Rscript")
        if not rscript:
            pytest.skip("Rscript not found on PATH")
        rcode = (
            "run_src <- readLines('run.R')\n"
            "run_src <- run_src[!grepl('^postde_main\\\\(\\\\)\\\\s*$', run_src)]\n"
            "srcfile <- tempfile(fileext = '.R')\n"
            "writeLines(run_src, srcfile)\n"
            "source(srcfile)\n"
            'source("common.R")\n'
            "config <- postde_config_defaults(list())\n"
            'cat("global", config$min_size, config$max_size, "\\n")\n'
            "f <- formals(postde_gsva_scores)\n"
            'cat("gsva", f$min_size, f$max_size, "\\n")\n'
            'source("enrichment.R")\n'
            "samples <- c('S1', 'S2', 'S3', 'S4')\n"
            "cond <- factor(c('A', 'A', 'B', 'B'), levels = c('A', 'B'))\n"
            "md <- data.frame(row.names = samples, condition = cond)\n"
            "design <- model.matrix(~condition, data = md)\n"
            "contrast <- c('(Intercept)' = 0, 'conditionB' = 1)\n"
            "entry <- list(\n"
            "    id = 'A-VS-B',\n"
            "    expression = matrix(rnorm(200), nrow = 50,\n"
            "        dimnames = list(paste0('G', 1:50), samples)),\n"
            "    design = design,\n"
            "    contrast = contrast\n"
            ")\n"
            "term2gene <- data.frame(term = c('T1', 'T1', 'T2', 'T2'),\n"
            "    gene = c('G1', 'G2', 'G3', 'G4'))\n"
            "config <- list(min_size = 10, max_size = 500,\n"
            "    gsva = list(enabled = TRUE))\n"
            "captured <- NULL\n"
            "postde_gsva_scores <- function(expr, gene_sets, method = 'gsva',\n"
            "    min_size = 10, max_size = 500) {\n"
            "    captured <<- c(min_size, max_size)\n"
            "    matrix(0, nrow = 1, ncol = ncol(expr),\n"
            "        dimnames = list('T1', colnames(expr)))\n"
            "}\n"
            "outdir <- tempfile()\n"
            "dir.create(outdir)\n"
            "invisible(postde_run_gsva(entry, term2gene, NULL, config, outdir))\n"
            'cat("effective", captured[1], captured[2], "\\n")\n'
        )
        rfile = tmp_path / "defaults.R"
        rfile.write_text(rcode)
        result = subprocess.run(
            [rscript, str(rfile)],
            cwd=str(REPO / "scripts" / "Analysis" / "PostDE"),
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stderr[-2000:]
        lines = [ln.split() for ln in result.stdout.strip().splitlines()]
        assert lines[0] == ["global", "10", "500"]
        assert lines[1] == ["gsva", "10", "500"]
        assert lines[2] == ["effective", "10", "500"]


class TestPostdeExecutable:
    @pytest.mark.parametrize("engine", ["smk", "nf"])
    @pytest.mark.parametrize("detool", ["deseq2", "edger"])
    def test_postde_runs_real_run_r(self, workdir, tmp_path, engine, detool):
        """Run the generated postde process/rule unmodified against the real
        run.R with an offline decoupleR ULM/MLM analysis on a tiny bundle."""
        rscript = shutil.which("Rscript")
        if not rscript:
            pytest.skip("Rscript not found on PATH")
        resdir = tmp_path / "res dir"
        resdir.mkdir()
        (resdir / "network.tsv").write_text(
            "source\ttarget\tmor\n"
            "TF1\tG1\t1\n"
            "TF1\tG2\t1\n"
            "TF2\tG3\t-1\n"
            "TF2\tG4\t-1\n"
        )
        analysis = {
            "padj": 0.05,
            "lfc": 1,
            "min_size": 2,
            "max_size": 500,
            "decoupler": {
                "enabled": True,
                "network": "res dir/network.tsv",
                "methods": ["ulm", "mlm"],
                "min_size": 2,
            },
        }
        cfg = _postde_cfg(tmp_path, analysis)
        subdir = _generate_de(engine, workdir, cfg)
        confo = json.load(
            open(workdir / subdir / f"allconditions_DE_{detool}_subconfig.json")
        )
        assert confo["POSTDE"]["enabled"] is True
        staged = confo["POSTDE"]["inputs"]
        launch = workdir / "launch dir"
        launch.mkdir()
        bundle_src = _write_real_bundle(rscript, launch, detool)
        if engine == "nf":
            subflow = (
                workdir / subdir / f"allconditions_DE_{detool}_subflow.nf"
            ).read_text()
            script = _extract_nf_script(subflow, "postde")
            script = script.replace("\\$", "$")
            script = script.replace("${params.gBINS}", str(REPO / "scripts"))
            shutil.copy2(bundle_src, launch / "bundle")
            os.symlink(workdir / staged, launch / "postde_inputs")
            result = subprocess.run(
                script, shell=True, cwd=str(launch), capture_output=True, text=True
            )
            assert result.returncode == 0, result.stderr[-2000:]
            outdir = launch / "postde"
            log = launch / "log"
        else:
            smk = (
                workdir / subdir / f"allconditions_DE_{detool}_subsnake.smk"
            ).read_text()
            body = _extract_smk_rule_body(smk, "postde")
            params_code, shell_tpl = _smk_params_shell(body)
            combo = f"{detool}_DE"
            ns = {"os": os, "combo": combo}
            exec(params_code, ns)
            fake_input = {
                "bundle": [str(bundle_src)],
                "staged": str(workdir / staged),
                "config": str(workdir / staged / "config.json"),
                "script": str(REPO / "scripts" / "Analysis" / "PostDE" / "run.R"),
                "common": str(REPO / "scripts" / "Analysis" / "PostDE" / "common.R"),
                "enrichment": str(
                    REPO / "scripts" / "Analysis" / "PostDE" / "enrichment.R"
                ),
                "regulatory": str(
                    REPO / "scripts" / "Analysis" / "PostDE" / "regulatory.R"
                ),
            }
            old = os.getcwd()
            os.chdir(launch)
            try:
                params = {
                    "bundle": ns["bundle"](
                        None, types.SimpleNamespace(**fake_input)
                    ),
                    "outdir": ns["outdir"](None),
                }
                log = os.path.abspath(
                    os.path.join("LOGS", combo, "DE", detool, "postde.log")
                )
            finally:
                os.chdir(old)
            shell = _format_shell(
                shell_tpl,
                input=types.SimpleNamespace(**fake_input),
                params=types.SimpleNamespace(**params),
                log=log,
            )
            result = subprocess.run(
                shell, shell=True, cwd=str(launch), capture_output=True, text=True
            )
            assert result.returncode == 0, result.stderr[-2000:]
            outdir = launch / "POSTDE" / combo
            log = launch / "LOGS" / combo / "DE" / detool / "postde.log"
        assert (outdir / "manifest.json").exists()
        assert (outdir / "report_data.rds").exists()
        entry = outdir / "A-VS-B"
        for name in (
            "decoupler_ulm_activities.tsv",
            "decoupler_ulm_differential.tsv",
            "decoupler_mlm_activities.tsv",
            "decoupler_mlm_differential.tsv",
            "decoupler_network_filtered.tsv",
            "decoupler_design.tsv",
            "decoupler_contrast.tsv",
        ):
            assert (entry / name).exists(), name
        assert log.exists()


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
