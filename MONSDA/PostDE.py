import hashlib
import json
import math
import os
import shutil

POSTDE_SECTION_KEYS = frozenset({"enabled", "config"})
POSTDE_TOP_KEYS = frozenset(
    {
        "padj",
        "lfc",
        "seed",
        "min_size",
        "max_size",
        "gprofiler",
        "clusterprofiler",
        "gsva",
        "decoupler",
        "dream",
        "plots",
        "report",
        "shiny",
    }
)
POSTDE_METHOD_KEYS = {
    "gprofiler": frozenset(
        {
            "enabled",
            "organism",
            "domain_scope",
            "allow_network",
            "sources",
            "correction_method",
            "background",
        }
    ),
    "clusterprofiler": frozenset(
        {"enabled", "term2gene", "term2name", "go_expand", "background"}
    ),
    "gsva": frozenset({"enabled", "method", "min_size", "max_size"}),
    "decoupler": frozenset(
        {
            "enabled",
            "network",
            "source_col",
            "target_col",
            "mor_col",
            "allow_network",
            "resource",
            "organism",
            "top",
            "min_size",
            "methods",
            "contrast_activity",
        }
    ),
    "dream": frozenset({"enabled", "metadata", "formula", "contrasts"}),
    "plots": frozenset({"enabled"}),
    "report": frozenset({"enabled"}),
    "shiny": frozenset({"enabled"}),
}
POSTDE_RESOURCE_KEYS = {
    "gprofiler": ("background",),
    "clusterprofiler": ("term2gene", "term2name", "background"),
    "decoupler": ("network",),
    "dream": ("metadata",),
}
POSTDE_BOOL_FLAGS = {
    "gprofiler": ("allow_network",),
    "clusterprofiler": ("go_expand",),
    "decoupler": ("allow_network", "contrast_activity"),
}
POSTDE_SUPPORTED_DE_TOOLS = frozenset({"deseq2", "edger"})
POSTDE_DOMAIN_SCOPES = ("annotated", "known", "custom", "custom_annotated")
POSTDE_CORRECTION_METHODS = (
    "g_SCS",
    "bonferroni",
    "fdr",
    "false_discovery_rate",
    "gSCS",
    "analytical",
)
POSTDE_GSVA_METHODS = ("gsva", "ssgsea")
POSTDE_DECOUPLER_RESOURCES = ("collectri", "progeny")
POSTDE_DECOUPLER_METHODS = ("ulm", "mlm")


def _is_number(value):
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(value)
    )


def _is_integral(value):
    return _is_number(value) and value == int(value)


def _validate_analysis_config(analysis, base_dir):
    unknown = set(analysis) - POSTDE_TOP_KEYS
    if unknown:
        raise ValueError(
            "POSTDE config has unknown top-level fields: "
            + ", ".join(sorted(unknown))
        )
    if "padj" in analysis and (
        not _is_number(analysis["padj"]) or not 0 < analysis["padj"] < 1
    ):
        raise ValueError("POSTDE config padj must be numeric in (0, 1)")
    if "lfc" in analysis and (
        not _is_number(analysis["lfc"]) or analysis["lfc"] < 0
    ):
        raise ValueError("POSTDE config lfc must be numeric >= 0")
    if "seed" in analysis and not _is_integral(analysis["seed"]):
        raise ValueError("POSTDE config seed must be a finite integer")
    if "min_size" in analysis and (
        not _is_integral(analysis["min_size"]) or analysis["min_size"] < 1
    ):
        raise ValueError("POSTDE config min_size must be an integer >= 1")
    if "max_size" in analysis and (
        not _is_integral(analysis["max_size"])
        or analysis["max_size"] < analysis.get("min_size", 1)
    ):
        raise ValueError("POSTDE config max_size must be an integer >= min_size")
    for method, keys in POSTDE_METHOD_KEYS.items():
        sub = analysis.get(method)
        if sub is None:
            continue
        if not isinstance(sub, dict):
            raise ValueError("POSTDE config " + method + " must be an object")
        unknown = set(sub) - keys
        if unknown:
            raise ValueError(
                "POSTDE config "
                + method
                + " has unknown fields: "
                + ", ".join(sorted(unknown))
            )
        if "enabled" in sub and not isinstance(sub["enabled"], bool):
            raise ValueError(
                "POSTDE config " + method + ".enabled must be a boolean"
            )
        for flag in POSTDE_BOOL_FLAGS.get(method, ()):
            if flag in sub and not isinstance(sub[flag], bool):
                raise ValueError(
                    "POSTDE config " + method + "." + flag + " must be a boolean"
                )
    for vis in ("plots", "report", "shiny"):
        if analysis.get(vis, {}).get("enabled"):
            raise ValueError(
                "POSTDE config "
                + vis
                + ".enabled=true is not available (reporting.R not yet provided)"
            )
    gp = analysis.get("gprofiler", {})
    if gp.get("enabled"):
        if not isinstance(gp.get("organism"), str) or not gp["organism"].strip():
            raise ValueError(
                "POSTDE config gprofiler.enabled=true requires gprofiler.organism"
            )
        if not gp.get("allow_network"):
            raise ValueError(
                "POSTDE config gprofiler.enabled=true requires gprofiler.allow_network=true"
            )
        if "domain_scope" in gp and gp["domain_scope"] not in POSTDE_DOMAIN_SCOPES:
            raise ValueError(
                "POSTDE config gprofiler.domain_scope must be one of "
                + ", ".join(POSTDE_DOMAIN_SCOPES)
            )
        if "correction_method" in gp and gp[
            "correction_method"
        ] not in POSTDE_CORRECTION_METHODS:
            raise ValueError(
                "POSTDE config gprofiler.correction_method must be one of "
                + ", ".join(POSTDE_CORRECTION_METHODS)
            )
    cp = analysis.get("clusterprofiler", {})
    gs = analysis.get("gsva", {})
    if (cp.get("enabled") or gs.get("enabled")) and not cp.get("term2gene"):
        raise ValueError(
            "POSTDE config clusterprofiler.term2gene required when clusterprofiler/gsva enabled"
        )
    if "method" in gs and gs["method"] not in POSTDE_GSVA_METHODS:
        raise ValueError(
            "POSTDE config gsva.method must be one of " + ", ".join(POSTDE_GSVA_METHODS)
        )
    if "min_size" in gs and (
        not _is_integral(gs["min_size"]) or gs["min_size"] < 1
    ):
        raise ValueError("POSTDE config gsva.min_size must be an integer >= 1")
    if "max_size" in gs and (
        not _is_integral(gs["max_size"])
        or gs["max_size"] < gs.get("min_size", analysis.get("min_size", 1))
    ):
        raise ValueError("POSTDE config gsva.max_size must be an integer >= min_size")
    dc = analysis.get("decoupler", {})
    if dc.get("enabled"):
        if not dc.get("network") and not dc.get("allow_network"):
            raise ValueError(
                "POSTDE config decoupler.enabled=true requires decoupler.network or decoupler.allow_network=true"
            )
        if "resource" in dc and dc["resource"] not in POSTDE_DECOUPLER_RESOURCES:
            raise ValueError(
                "POSTDE config decoupler.resource must be one of "
                + ", ".join(POSTDE_DECOUPLER_RESOURCES)
            )
        if "methods" in dc:
            methods = (
                dc["methods"] if isinstance(dc["methods"], list) else [dc["methods"]]
            )
            if not methods or any(
                m not in POSTDE_DECOUPLER_METHODS for m in methods
            ):
                raise ValueError(
                    "POSTDE config decoupler.methods must be a list of "
                    + ", ".join(POSTDE_DECOUPLER_METHODS)
                )
        if "min_size" in dc and (
            not _is_integral(dc["min_size"]) or dc["min_size"] < 1
        ):
            raise ValueError("POSTDE config decoupler.min_size must be an integer >= 1")
        if "top" in dc and (not _is_integral(dc["top"]) or dc["top"] < 1):
            raise ValueError("POSTDE config decoupler.top must be an integer >= 1")
    dm = analysis.get("dream", {})
    if dm.get("enabled"):
        if not isinstance(dm.get("metadata"), str) or not dm["metadata"].strip():
            raise ValueError(
                "POSTDE config dream.enabled=true requires dream.metadata"
            )
        if not isinstance(dm.get("formula"), str) or not dm["formula"].strip():
            raise ValueError("POSTDE config dream.enabled=true requires dream.formula")
    for method, keys in POSTDE_RESOURCE_KEYS.items():
        sub = analysis.get(method, {})
        for key in keys:
            path = sub.get(key)
            if path is None:
                continue
            if not isinstance(path, str) or not path.strip():
                raise ValueError(
                    "POSTDE config " + method + "." + key + " must be a non-empty path"
                )
            full = path if os.path.isabs(path) else os.path.join(base_dir, path)
            if not os.path.isfile(full):
                raise ValueError(
                    "POSTDE config resource file not found: "
                    + method
                    + "."
                    + key
                    + " -> "
                    + full
                )


def load_postde_config(section, config=None):
    if section is None:
        return None
    if not isinstance(section, dict):
        raise ValueError("POSTDE section must be an object")
    unknown = set(section) - POSTDE_SECTION_KEYS
    if unknown:
        raise ValueError(
            "POSTDE section has unknown fields: " + ", ".join(sorted(unknown))
        )
    if "enabled" not in section:
        raise ValueError("POSTDE section missing required field 'enabled'")
    if not isinstance(section["enabled"], bool):
        raise ValueError("POSTDE.enabled must be a boolean")
    if not section["enabled"]:
        return {"enabled": False, "config": section.get("config", "")}
    if not isinstance(section.get("config"), str) or not section["config"].strip():
        raise ValueError("POSTDE.enabled=true requires a non-empty POSTDE.config path")
    config_path = os.path.abspath(section["config"])
    if not os.path.isfile(config_path):
        raise ValueError("POSTDE.config file not found: " + config_path)
    try:
        with open(config_path) as fh:
            analysis = json.load(fh)
    except (OSError, ValueError) as err:
        raise ValueError("POSTDE.config is not readable JSON: " + str(err))
    if not isinstance(analysis, dict):
        raise ValueError("POSTDE.config must contain a JSON object")
    _validate_analysis_config(analysis, os.path.dirname(config_path))
    if config is not None:
        workflows = [
            w.strip()
            for w in str(config.get("WORKFLOWS", "")).split(",")
            if w.strip()
        ]
        if "DE" not in workflows:
            raise ValueError("POSTDE.enabled=true requires DE in WORKFLOWS")
        de_tools = config.get("DE", {}).get("TOOLS", {})
        if isinstance(de_tools, dict):
            unsupported = set(de_tools) - POSTDE_SUPPORTED_DE_TOOLS
            if unsupported:
                raise ValueError(
                    "POSTDE.enabled=true does not support DE tools: "
                    + ", ".join(sorted(unsupported))
                )
    return {"enabled": True, "config": config_path}


def _collect_resources(analysis, base_dir):
    resources = []
    for method, keys in POSTDE_RESOURCE_KEYS.items():
        sub = analysis.get(method, {})
        for key in keys:
            path = sub.get(key)
            if path is None:
                continue
            full = path if os.path.isabs(path) else os.path.join(base_dir, path)
            resources.append((method, key, full))
    return resources


def _staged_name(method, key, full):
    ext = os.path.splitext(os.path.basename(full))[1]
    return method + "_" + key + ext


def _same_content(path_a, path_b):
    """Compare file contents directly (readbytes), immune to stat-based
    caching that can miss same-size/same-mtime content mutations."""
    with open(path_a, "rb") as fh_a, open(path_b, "rb") as fh_b:
        return fh_a.read() == fh_b.read()


def _same_text(path, text):
    """Compare a file's bytes against an expected text string."""
    with open(path, "rb") as fh:
        return fh.read() == text.encode("utf-8")


def prepare_postde(config, subdir):
    section = load_postde_config(config.get("POSTDE"), config=config)
    if section is None or not section.get("enabled"):
        return None
    config_path = section["config"]
    with open(config_path) as fh:
        analysis = json.load(fh)
    base_dir = os.path.dirname(config_path)
    resources = _collect_resources(analysis, base_dir)
    digest = hashlib.sha256()
    digest.update(json.dumps(analysis, sort_keys=True).encode("utf-8"))
    for method, key, full in sorted(resources):
        with open(full, "rb") as fh:
            digest.update(fh.read())
    staged = os.path.join(os.path.abspath(subdir), "POSTDE_" + digest.hexdigest()[:12])
    names = {
        (method, key): _staged_name(method, key, full)
        for method, key, full in resources
    }
    staged_analysis = json.loads(json.dumps(analysis))
    for method, key, full in resources:
        staged_analysis[method][key] = names[(method, key)]
    expected_config = json.dumps(staged_analysis, indent=2)
    if os.path.isdir(staged):
        staged_config = os.path.join(staged, "config.json")
        config_ok = os.path.isfile(staged_config) and _same_text(
            staged_config, expected_config
        )
        resources_ok = all(
            os.path.isfile(os.path.join(staged, names[(m, k)]))
            and _same_content(full, os.path.join(staged, names[(m, k)]))
            for m, k, full in resources
        )
        if config_ok and resources_ok:
            return staged
    os.makedirs(staged, exist_ok=True)
    for method, key, full in resources:
        dst = os.path.join(staged, names[(method, key)])
        if not os.path.exists(dst) or not _same_content(full, dst):
            shutil.copy2(full, dst)
    staged_config = os.path.join(staged, "config.json")
    if not os.path.isfile(staged_config) or not _same_text(
        staged_config, expected_config
    ):
        with open(staged_config, "w") as fh:
            fh.write(expected_config)
    return staged
