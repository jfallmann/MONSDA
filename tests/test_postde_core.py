import json
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RSCRIPT = shutil.which("Rscript")
pytestmark = pytest.mark.skipif(
    RSCRIPT is None, reason="Rscript is not on PATH; activate the postde Conda environment"
)
POSTDE_DIR = ROOT / "scripts" / "Analysis" / "PostDE"
FIXTURE = ROOT / "tests" / "postde_core.R"

REQUIRED_SUMMARY_CHECKS = [
    "deseq2_rownames_aligned",
    "deseq2_gene_id_matches",
    "edger_rownames_aligned",
    "edger_gene_id_matches",
    "edger_stat_finite",
    "edger_stat_na_noF",
    "capture_rejects_misaligned",
    "capture_counts_superset_ok",
    "diff_feature_names_preserved",
    "decoupler_network_meta_exists",
    "decoupler_filtered_tsv_exists",
    "decoupler_gated_without_network",
    "mlm_rank_deficiency_rejected",
    "gomap_identical",
    "gomap_unknown_term_rejected",
    "term2gene_strict",
    "term2name_strict",
    "term2gene_extra_col_rejected",
    "ora_separate_files",
    "enrichment_provenance_exists",
    "gsea_ranked_saved",
    "dream_varpart_status_file",
    "dream_varpart_ok_path",
    "dream_safe_formula_ok",
    "dream_unsafe_formula_rejected",
    "dream_safe_contrast_ok",
    "dream_unsafe_contrast_rejected",
    "empty_tsv_header_first",
    "empty_tsv_status_sidecar",
    "sig_threshold_strict",
    "sig_never_select_na_effect",
    "singleton_no_residual_df_rejected",
    "design_named_contrast_aligned",
    "design_zero_intercept_rejected",
    "sanitize_dot_rejected",
    "sanitize_collision_detected",
]


def run_rscript(args, cwd=None, check=True, timeout=1200):
    cmd = [RSCRIPT] + [str(a) for a in args]
    res = subprocess.run(cmd, cwd=cwd or ROOT, capture_output=True, text=True, timeout=timeout)
    if check and res.returncode != 0:
        raise AssertionError(
            "Rscript failed with exit %d\nstdout:\n%s\nstderr:\n%s"
            % (res.returncode, res.stdout[-3000:], res.stderr[-3000:])
        )
    return res


@pytest.fixture(scope="module")
def fixture_dir(tmp_path_factory):
    out = tmp_path_factory.mktemp("postde_fixture")
    res = run_rscript([FIXTURE, out])
    assert "POSTDE_FIXTURE_OK" in res.stdout
    return out


@pytest.fixture(scope="module")
def summary(fixture_dir):
    with open(fixture_dir / "summary.json") as fh:
        return json.load(fh)


def test_fixture_runs_and_all_checks_pass(summary):
    for key in REQUIRED_SUMMARY_CHECKS:
        assert summary.get(key) is True, "summary check %s = %r" % (key, summary.get(key))
    assert summary["diff_columns"] == "feature_id,score_difference,pvalue,padj,stat,mean_score,SE"
    assert summary["gsea_rank_unavailable_status"] == "rank unavailable"
    assert summary["dream_varpart_status"] == "error"
    assert summary["dream_skipped_normalized"].startswith("skipped")
    assert summary["gsva_max_diff"] < 1e-6
    assert summary["diff_score_difference_max_diff"] < 1e-6
    assert summary["diff_stat_max_diff"] < 1e-6
    assert summary["ulm_max_diff"] < 1e-6
    assert summary["mlm_max_diff"] < 1e-6
    assert summary["dream_max_logfc_diff"] < 0.1
    assert summary["gomap_expanded_rows"] > 4


def test_bundles_written(fixture_dir):
    assert (fixture_dir / "DE_deseq2_test_postde.rds").exists()
    assert (fixture_dir / "DE_edger_test_postde.rds").exists()


def test_empty_tsv_header_first(fixture_dir):
    lines = (fixture_dir / "empty_result.tsv").read_text().splitlines()
    assert lines[0] == "gene_id\tlogFC\tpvalue\tpadj\tstat\tmean"
    status = json.loads((fixture_dir / "empty_result.tsv.status.json").read_text())
    assert status["status"] == "no significant genes"


def _write_config(path, data):
    path.write_text(json.dumps(data))
    return path


def test_run_cli_end_to_end(fixture_dir, tmp_path):
    config = {
        "padj": 0.99,
        "lfc": 0,
        "min_size": 1,
        "max_size": 500,
        "seed": 1,
        "clusterprofiler": {
            "enabled": True,
            "term2gene": str(fixture_dir / "term2gene.tsv"),
            "term2name": str(fixture_dir / "term2name.tsv"),
        },
        "decoupler": {
            "enabled": True,
            "network": str(fixture_dir / "test_network.tsv"),
            "min_size": 1,
            "methods": ["ulm", "mlm"],
        },
    }
    config_path = _write_config(tmp_path / "config.json", config)
    out = tmp_path / "out"
    res = run_rscript(
        [
            POSTDE_DIR / "run.R",
            "--bundle",
            fixture_dir / "DE_deseq2_test_postde.rds",
            "--config",
            config_path,
            "--output",
            out,
        ]
    )
    assert res.returncode == 0
    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["output"] == "."
    assert manifest["bundle_checksum"]
    assert manifest["config_checksum"]
    assert manifest["sessionInfo"]
    entry = manifest["entries"]["A_vs_B"]
    assert entry["dir"] == "A_vs_B"
    assert entry["analyses"]["enrichment"]["enricher_all"]["result"] == "A_vs_B/enricher_all_result.tsv"
    assert entry["analyses"]["regulatory"]["decoupler"]["ulm"]["activities"] == "A_vs_B/decoupler_ulm_activities.tsv"
    assert (out / "A_vs_B" / "enricher_all_result.tsv").exists()
    assert (out / "A_vs_B" / "decoupler_ulm_activities.tsv").exists()
    assert (out / "A_vs_B" / "decoupler_mlm_activities.tsv").exists()
    assert (out / "A_vs_B" / "gsea_result.tsv").exists()
    assert (out / "report_data.rds").exists()


def test_gprofiler_requires_allow_network(fixture_dir, tmp_path):
    config = _write_config(
        tmp_path / "config_gp.json",
        {"gprofiler": {"enabled": True, "organism": "hsapiens"}},
    )
    res = run_rscript(
        [
            POSTDE_DIR / "run.R",
            "--bundle",
            fixture_dir / "DE_deseq2_test_postde.rds",
            "--config",
            config,
            "--output",
            tmp_path / "out_gp",
        ],
        check=False,
    )
    assert res.returncode != 0
    assert "allow_network" in res.stderr


def test_gprofiler_mock_injection(fixture_dir, tmp_path):
    script = tmp_path / "mock_gprofiler.R"
    template = """
source("__COMMON__")
source("__ENRICH__")
mock_gost <- function(query, organism, domain_scope, sources, significant,
                      user_threshold, correction_method, custom_bg = NULL) {
    list(
        result = data.frame(
            query = "q1", significant = TRUE, p_value = 0.01,
            term_size = 10, query_size = 5, intersection_size = 3,
            effective_domain_size = 100, source = "GO:BP",
            term_id = "GO:0006915", term_name = "apoptotic process",
            parents = I(list(c("GO:0008150"))), stringsAsFactors = FALSE
        ),
        meta = list(version = "mock", organism = organism)
    )
}
entry <- list(results = data.frame(
    gene_id = paste0("g", 1:20), logFC = rnorm(20), pvalue = runif(20),
    padj = runif(20), stat = rnorm(20), stringsAsFactors = FALSE))
config <- list(padj = 0.5, lfc = 0,
    gprofiler = list(enabled = TRUE, organism = "hsapiens",
                     allow_network = TRUE, domain_scope = "custom",
                     sources = list("GO:BP"), correction_method = "g_SCS"))
out <- postde_run_gprofiler(entry, config, "__OUT__", gost_fun = mock_gost)
stopifnot(out$all$status == "ok")
stopifnot(file.exists(out$all$result))
stopifnot(file.exists(out$all$meta))
stopifnot(file.exists(out$all$raw))
res <- read.delim(out$all$result, check.names = FALSE)
stopifnot("term_id" %in% colnames(res))
cat("MOCK_GPROFILER_OK\n")
"""
    script.write_text(
        template.replace("__COMMON__", str(POSTDE_DIR / "common.R"))
        .replace("__ENRICH__", str(POSTDE_DIR / "enrichment.R"))
        .replace("__OUT__", str(tmp_path))
    )
    res = run_rscript([script])
    assert "MOCK_GPROFILER_OK" in res.stdout


def test_config_validation_rejects_bad_padj(fixture_dir, tmp_path):
    config = _write_config(tmp_path / "config_bad.json", {"padj": 1.5})
    res = run_rscript(
        [
            POSTDE_DIR / "run.R",
            "--bundle",
            fixture_dir / "DE_deseq2_test_postde.rds",
            "--config",
            config,
            "--output",
            tmp_path / "out_bad",
        ],
        check=False,
    )
    assert res.returncode != 0
    assert "padj" in res.stderr


def test_empty_bundle_rejected(tmp_path):
    bundle = tmp_path / "empty_bundle.rds"
    res = run_rscript(
        [
            "-e",
            "saveRDS(list(schema_version = 1L, engine = 'test', contrasts = list()), '%s')" % bundle,
        ]
    )
    assert res.returncode == 0
    config = _write_config(tmp_path / "config_empty.json", {})
    res = run_rscript(
        [
            POSTDE_DIR / "run.R",
            "--bundle",
            bundle,
            "--config",
            config,
            "--output",
            tmp_path / "out_empty",
        ],
        check=False,
    )
    assert res.returncode != 0
    assert "no contrast entries" in res.stderr


def test_duplicate_sanitized_ids_rejected(tmp_path):
    bundle = tmp_path / "dup_bundle.rds"
    script = tmp_path / "make_dup_bundle.R"
    entry_r = (
        "list(id = 'a/b', A = 'A', B = 'B', normalized = FALSE, metadata = list(), "
        "counts = list(), expression = list(), design = list(), contrast = list(), "
        "results = list(), mean_scale = 'baseMean')"
    )
    script.write_text(
        "saveRDS(list(schema_version = 1L, engine = 'test', "
        "contrasts = list(e1 = %s, e2 = %s)), '%s')\n" % (entry_r, entry_r, bundle)
    )
    res = run_rscript([script])
    assert res.returncode == 0
    config = _write_config(tmp_path / "config_dup.json", {})
    res = run_rscript(
        [
            POSTDE_DIR / "run.R",
            "--bundle",
            bundle,
            "--config",
            config,
            "--output",
            tmp_path / "out_dup",
        ],
        check=False,
    )
    assert res.returncode != 0
    assert "duplicate" in res.stderr.lower() or "unsafe" in res.stderr.lower()
