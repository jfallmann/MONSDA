import copy
import json
import os
import sys
from typing import Any, Dict, List, Optional

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse
from pydantic import BaseModel, Field
from snakemake.common.configfile import load_configfile

from .Params import samplesheet_to_settings

TEMPLATE_FILE = "template_base_commented.json"
NONE_WORKFLOW_KEYS = ["WORKFLOWS", "BINS", "MAXTHREADS", "SETTINGS", "VERSION"]

app = FastAPI(title="MONSDA Configurator Web")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


class SamplesheetRequest(BaseModel):
    samplesheet_path: str


class BuildConfigRequest(BaseModel):
    config_name: str
    output_dir: str
    workflows: List[str] = Field(default_factory=list)
    tools: Dict[str, List[str]] = Field(default_factory=dict)
    maxthreads: str = "16"
    settings: Optional[Dict[str, Any]] = None
    samplesheet_path: Optional[str] = None


class ConditionFiles(BaseModel):
    """Per-condition file selections from the wizard."""

    condition: str  # slash-separated condition path e.g. "Ecoli/WT"
    fastq_dir: str = ""  # directory containing FASTQ files for this condition
    fastq_files: List[str] = Field(default_factory=list)  # explicit file paths
    sequencing: str = "paired"  # paired | single
    reference: str = ""  # absolute path to genome FASTA
    gtf: str = ""  # absolute path to GTF annotation
    gff: str = ""  # absolute path to GFF annotation (optional)
    decoy: str = ""  # absolute path to decoy file (optional)


class ProjectRequest(BaseModel):
    project_dir: str
    project_name: str = "monsda"
    condition_files: List[ConditionFiles] = Field(default_factory=list)
    config: Optional[Dict[str, Any]] = None
    settings: Optional[Dict[str, Any]] = None


class SaveConfigRequest(BaseModel):
    config_name: str
    output_dir: str
    config: Dict[str, Any]

def _template_candidates() -> List[str]:
    base_dir = os.path.abspath(os.path.dirname(__file__))
    pythonversion = f"python{sys.version_info.major}.{sys.version_info.minor}"
    install_share = base_dir.replace(
        os.sep.join(["lib", pythonversion, "site-packages", "MONSDA"]), "share"
    )

    candidates = [
        # explicit override if user/admin wants to force a template location
        os.environ.get("MONSDA_TEMPLATE_PATH", ""),
        # local checkout layout (repo root/configs)
        os.path.abspath(os.path.join(base_dir, "..", "configs", TEMPLATE_FILE)),
        # package-relative layout if configs were bundled next to package
        os.path.abspath(os.path.join(base_dir, "configs", TEMPLATE_FILE)),
        # Configurator.py install layout logic: <prefix>/share/MONSDA/configs
        os.path.abspath(
            os.path.join(install_share, "MONSDA", "configs", TEMPLATE_FILE)
        ),
        # generic venv/conda prefix share path
        os.path.abspath(
            os.path.join(sys.prefix, "share", "MONSDA", "configs", TEMPLATE_FILE)
        ),
        # last-resort: run directory configs
        os.path.abspath(os.path.join(os.getcwd(), "configs", TEMPLATE_FILE)),
    ]

    # de-duplicate while preserving order and dropping empty values
    deduped: List[str] = []
    for c in candidates:
        c = (c or "").strip()
        if c and c not in deduped:
            deduped.append(c)
    return deduped


def resolve_template_path() -> str:
    for path in _template_candidates():
        if os.path.isfile(path):
            return path
    tried = "\n - ".join(_template_candidates())
    raise FileNotFoundError(
        "Could not find template_base_commented.json. Tried:\n - " + tried
    )


def load_template() -> Dict[str, Any]:
    return load_configfile(resolve_template_path())


def strip_comments(d: Any) -> Any:
    if isinstance(d, dict):
        return {k: strip_comments(v) for k, v in d.items() if k != "comment"}
    if isinstance(d, list):
        return [strip_comments(x) for x in d]
    return d


def set_in_path(root: Dict[str, Any], path: List[str], value: Any) -> None:
    node = root
    for key in path[:-1]:
        if key not in node or not isinstance(node[key], dict):
            node[key] = {}
        node = node[key]
    node[path[-1]] = value


def get_condition_paths_from_settings(settings: Dict[str, Any]) -> List[List[str]]:
    out: List[List[str]] = []

    def _walk(node: Any, path: List[str]) -> None:
        if not isinstance(node, dict):
            return
        if "SAMPLES" in node and isinstance(node.get("SAMPLES"), list):
            out.append(path.copy())
            return
        for k, v in node.items():
            if isinstance(v, dict):
                _walk(v, path + [k])

    _walk(settings, [])
    return out


def build_workflow_block(
    workflow_name: str,
    workflow_template: Dict[str, Any],
    condition_paths: List[List[str]],
    selected_tools: List[str],
) -> Dict[str, Any]:
    block: Dict[str, Any] = {}

    all_tools = workflow_template.get("TOOLS", {})
    if not selected_tools:
        selected_tools = list(all_tools.keys())

    if all_tools:
        block["TOOLS"] = {k: all_tools[k] for k in selected_tools if k in all_tools}

    # carry workflow-level defaults if present
    for passthrough in ["FEATURES", "CUTOFFS", "COMPARABLE", "EXCLUDE"]:
        if passthrough in workflow_template:
            block[passthrough] = copy.deepcopy(workflow_template[passthrough])

    # populate per-condition tool settings like CLI configurator does
    for cond_path in condition_paths:
        for tool in selected_tools:
            tool_def = workflow_template.get(tool)
            if not isinstance(tool_def, dict):
                continue
            tool_def_no_comment = strip_comments(tool_def)
            if not tool_def_no_comment:
                continue
            set_in_path(block, cond_path + [tool], copy.deepcopy(tool_def_no_comment))

    return block


def build_config(req: BuildConfigRequest) -> Dict[str, Any]:
    template = strip_comments(load_template())

    if req.settings is not None:
        settings = req.settings
    elif req.samplesheet_path:
        samplesheet_path = os.path.abspath(req.samplesheet_path)
        if not os.path.isfile(samplesheet_path):
            raise HTTPException(
                status_code=400,
                detail=f"Samplesheet does not exist: {samplesheet_path}",
            )
        settings = samplesheet_to_settings(samplesheet_path)
    else:
        raise HTTPException(
            status_code=400,
            detail="Provide either settings JSON or samplesheet_path.",
        )

    condition_paths = get_condition_paths_from_settings(settings)
    if not condition_paths:
        raise HTTPException(
            status_code=400,
            detail="No condition leaves with SAMPLES found in SETTINGS.",
        )

    workflows = req.workflows or []
    if not workflows:
        raise HTTPException(
            status_code=400,
            detail="At least one workflow must be selected.",
        )

    invalid = [w for w in workflows if w not in template or w in NONE_WORKFLOW_KEYS]
    if invalid:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown workflow(s): {', '.join(invalid)}",
        )

    final_config: Dict[str, Any] = {
        "WORKFLOWS": ", ".join(workflows),
        "BINS": template.get("BINS", ""),
        "MAXTHREADS": str(req.maxthreads),
        "VERSION": template.get("VERSION", ""),
        "SETTINGS": settings,
    }

    for wf in workflows:
        wf_template = template.get(wf, {})
        final_config[wf] = build_workflow_block(
            wf,
            wf_template,
            condition_paths,
            req.tools.get(wf, []),
        )

    return final_config


@app.get("/template", response_model=Dict[str, Any])
def get_template() -> Dict[str, Any]:
    try:
        return load_template()
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to load template: {e}")


@app.get("/template/fields", response_model=Dict[str, Any])
def get_template_fields() -> Dict[str, Any]:
    try:
        return strip_comments(load_template())
    except Exception as e:
        raise HTTPException(
            status_code=500, detail=f"Failed to load template fields: {e}"
        )


def _normalize_path(path: str) -> str:
    if path:
        return os.path.abspath(os.path.expanduser(path.strip()))
    return os.getcwd()


@app.get("/fs/roots", response_model=Dict[str, Any])
def fs_roots() -> Dict[str, Any]:
    cwd = os.getcwd()
    home = os.path.expanduser("~")
    roots = []
    for p in [cwd, home, os.path.sep]:
        if p and p not in roots and os.path.isdir(p):
            roots.append(p)
    return {"roots": roots}


@app.get("/fs/list", response_model=Dict[str, Any])
def fs_list(path: str = "", mode: str = "dirs") -> Dict[str, Any]:
    mode = (mode or "dirs").lower()
    if mode not in {"dirs", "all", "samplesheets"}:
        raise HTTPException(
            status_code=400, detail="mode must be dirs, all, or samplesheets"
        )

    current = _normalize_path(path)
    if not os.path.isdir(current):
        raise HTTPException(status_code=400, detail=f"Not a directory: {current}")

    try:
        entries = list(os.scandir(current))
    except PermissionError:
        raise HTTPException(status_code=403, detail=f"Permission denied: {current}")

    dirs = []
    files = []
    for e in sorted(entries, key=lambda x: x.name.lower()):
        if e.name in {".", ".."}:
            continue
        if e.is_dir(follow_symlinks=True):
            dirs.append({"name": e.name, "path": e.path})
            continue
        if not e.is_file(follow_symlinks=True):
            continue
        if mode in {"all", "samplesheets"}:
            files.append({"name": e.name, "path": e.path})

    parent = (
        os.path.dirname(current.rstrip(os.sep))
        if current != os.path.sep
        else os.path.sep
    )
    return {
        "current": current,
        "parent": parent,
        "dirs": dirs,
        "files": files,
        "mode": mode,
    }


@app.post("/samplesheet/parse", response_model=Dict[str, Any])
def parse_samplesheet(req: SamplesheetRequest) -> Dict[str, Any]:
    path = os.path.abspath(req.samplesheet_path.strip())
    if not os.path.isfile(path):
        raise HTTPException(status_code=400, detail=f"Samplesheet not found: {path}")
    settings = samplesheet_to_settings(path)
    return {
        "samplesheet": path,
        "settings": settings,
        "conditions": [
            "/".join(p) for p in get_condition_paths_from_settings(settings)
        ],
    }


@app.post("/config/preview", response_model=Dict[str, Any])
def preview_config(req: BuildConfigRequest) -> Dict[str, Any]:
    return {"config": build_config(req)}


@app.post("/config/save", response_model=Dict[str, Any])
def save_config(req: SaveConfigRequest) -> Dict[str, Any]:
    config_name = req.config_name.strip()
    if not config_name or any(c in config_name for c in "/\\"):
        raise HTTPException(status_code=400, detail="Invalid config name.")

    output_dir = os.path.abspath(req.output_dir.strip())
    if not os.path.isdir(output_dir):
        raise HTTPException(status_code=400, detail="Output directory does not exist.")

    path = os.path.join(output_dir, f"config_{config_name}.json")
    with open(path, "w") as fh:
        json.dump(req.config, fh, indent=4)

    return {"message": "Config written.", "path": path}


@app.post("/project/create", response_model=Dict[str, Any])
def create_project(req: ProjectRequest) -> Dict[str, Any]:
    base_dir = os.path.abspath(req.project_dir.strip())
    project_name = req.project_name.strip() or "monsda"

    # Project name becomes a subdirectory under the chosen path
    project_dir = os.path.join(base_dir, project_name)
    os.makedirs(project_dir, exist_ok=True)

    # FASTQ directly under project dir (no extra project_name level)
    fastq_dir = os.path.join(project_dir, "FASTQ")
    gen_dir = os.path.join(project_dir, "GENOMES")
    os.makedirs(fastq_dir, exist_ok=True)
    os.makedirs(gen_dir, exist_ok=True)

    settings = req.settings or (req.config.get("SETTINGS") if req.config else None)
    linked_files: List[str] = []
    warnings: List[str] = []

    # Build lookup: condition_path -> ConditionFiles
    cond_lookup: Dict[str, ConditionFiles] = {}
    for cf in req.condition_files:
        cond_lookup[cf.condition] = cf

    if settings and isinstance(settings, dict):
        condition_paths = get_condition_paths_from_settings(settings)
        for cond_path in condition_paths:
            cond_key = "/".join(cond_path)
            cond_info = cond_lookup.get(cond_key)

            # Create condition subdirectories under FASTQ
            cond_dir = os.path.join(fastq_dir, *cond_path)
            os.makedirs(cond_dir, exist_ok=True)

            # Walk to the leaf node to find SAMPLES
            node = settings
            for key in cond_path:
                node = node.get(key, {})

            samples = node.get("SAMPLES", [])

            # Link FASTQ files: prefer explicit file list, fall back to directory+sample matching
            if cond_info and cond_info.fastq_files:
                # Explicit file selection mode - link exactly the files the user chose
                for fpath in cond_info.fastq_files:
                    src_file = os.path.abspath(fpath)
                    if not os.path.isfile(src_file):
                        warnings.append(f"FASTQ file not found: {fpath}")
                        continue
                    dst = os.path.join(cond_dir, os.path.basename(src_file))
                    if not os.path.exists(dst):
                        os.symlink(os.path.realpath(src_file), dst)
                        linked_files.append(dst)
                    # Handle paired-end mate
                    if cond_info.sequencing == "paired":
                        basename = os.path.basename(src_file)
                        if "_R1" in basename:
                            mate_name = basename.replace("_R1", "_R2")
                        elif "_1." in basename:
                            mate_name = basename.replace("_1.", "_2.")
                        else:
                            mate_name = None
                        if mate_name:
                            mate_src = os.path.join(
                                os.path.dirname(src_file), mate_name
                            )
                            mate_dst = os.path.join(cond_dir, mate_name)
                            if os.path.isfile(mate_src) and not os.path.exists(
                                mate_dst
                            ):
                                os.symlink(os.path.realpath(mate_src), mate_dst)
                                linked_files.append(mate_dst)

            elif cond_info and cond_info.fastq_dir and isinstance(samples, list):
                # Directory mode - find files matching sample names from samplesheet
                src_dir = os.path.abspath(cond_info.fastq_dir)
                if os.path.isdir(src_dir):
                    for sample in samples:
                        # Find matching files (sample name may or may not have extension)
                        matched = _find_fastq_files(src_dir, sample)
                        if not matched:
                            warnings.append(
                                f"No FASTQ file found for sample '{sample}' in {src_dir}"
                            )
                            continue
                        for src_file in matched:
                            dst = os.path.join(cond_dir, os.path.basename(src_file))
                            if not os.path.exists(dst):
                                os.symlink(os.path.realpath(src_file), dst)
                                linked_files.append(dst)
                            # Handle paired-end mate
                            if cond_info.sequencing == "paired":
                                basename = os.path.basename(src_file)
                                if "_R1" in basename:
                                    mate_name = basename.replace("_R1", "_R2")
                                elif "_1." in basename:
                                    mate_name = basename.replace("_1.", "_2.")
                                elif "_R2" in basename:
                                    mate_name = basename.replace("_R2", "_R1")
                                elif "_2." in basename:
                                    mate_name = basename.replace("_2.", "_1.")
                                else:
                                    mate_name = None
                                if mate_name:
                                    mate_src = os.path.join(src_dir, mate_name)
                                    mate_dst = os.path.join(cond_dir, mate_name)
                                    if os.path.isfile(mate_src) and not os.path.exists(
                                        mate_dst
                                    ):
                                        os.symlink(os.path.realpath(mate_src), mate_dst)
                                        linked_files.append(mate_dst)
                else:
                    warnings.append(
                        f"FASTQ directory not found for condition '{cond_key}': {src_dir}"
                    )

            # Link genome files: REFERENCE, DECOY, GTF, GFF into GENOMES/
            ref_path = (cond_info.reference if cond_info else "") or node.get(
                "REFERENCE", ""
            )
            if ref_path and os.path.isfile(ref_path):
                dst = os.path.join(gen_dir, os.path.basename(ref_path))
                if not os.path.exists(dst):
                    os.symlink(os.path.realpath(ref_path), dst)
                    linked_files.append(dst)
                # Update settings path to relative
                rel = os.path.relpath(dst, start=project_dir)
                node["REFERENCE"] = rel
            elif ref_path:
                warnings.append(f"REFERENCE not found: {ref_path}")

            decoy_path = (cond_info.decoy if cond_info else "") or node.get("DECOY", "")
            if decoy_path and os.path.isfile(decoy_path):
                dst = os.path.join(gen_dir, os.path.basename(decoy_path))
                if not os.path.exists(dst):
                    os.symlink(os.path.realpath(decoy_path), dst)
                    linked_files.append(dst)
                rel = os.path.relpath(dst, start=project_dir)
                node["DECOY"] = rel

            gtf_path = (cond_info.gtf if cond_info else "") or ""
            anno = node.get("ANNOTATION", {})
            if isinstance(anno, dict):
                gtf_path = gtf_path or anno.get("GTF", "")
            if gtf_path and os.path.isfile(gtf_path):
                dst = os.path.join(gen_dir, os.path.basename(gtf_path))
                if not os.path.exists(dst):
                    os.symlink(os.path.realpath(gtf_path), dst)
                    linked_files.append(dst)
                rel = os.path.relpath(dst, start=project_dir)
                if isinstance(anno, dict):
                    anno["GTF"] = rel
                    node["ANNOTATION"] = anno
            elif gtf_path:
                warnings.append(f"GTF not found: {gtf_path}")

            gff_path = (cond_info.gff if cond_info else "") or ""
            if isinstance(anno, dict):
                gff_path = gff_path or anno.get("GFF", "")
            if gff_path and os.path.isfile(gff_path):
                dst = os.path.join(gen_dir, os.path.basename(gff_path))
                if not os.path.exists(dst):
                    os.symlink(os.path.realpath(gff_path), dst)
                    linked_files.append(dst)
                rel = os.path.relpath(dst, start=project_dir)
                if isinstance(anno, dict):
                    anno["GFF"] = rel
                    node["ANNOTATION"] = anno
            elif gff_path:
                warnings.append(f"GFF not found: {gff_path}")

    # Write config if provided (with updated relative paths in SETTINGS)
    config_path = ""
    if req.config:
        if settings:
            req.config["SETTINGS"] = settings
        config_file = f"config_{project_name}.json"
        config_path = os.path.join(project_dir, config_file)
        with open(config_path, "w") as fh:
            json.dump(req.config, fh, indent=4)

    return {
        "message": "Project created.",
        "path": project_dir,
        "config_path": config_path,
        "linked_files": len(linked_files),
        "warnings": warnings,
    }


def _find_fastq_files(src_dir: str, sample_name: str) -> List[str]:
    """Find FASTQ files in src_dir matching a sample name.

    Matches: exact filename, or sample_name*.fastq.gz / .fq.gz / .fastq / .fq
    Only returns R1 (or unpaired) to avoid double-linking.
    """
    import glob as _glob

    candidates: List[str] = []
    # Try exact match first
    exact = os.path.join(src_dir, sample_name)
    if os.path.isfile(exact):
        return [exact]

    # Try common FASTQ extensions
    for ext in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
        path = os.path.join(src_dir, sample_name + ext)
        if os.path.isfile(path):
            candidates.append(path)
        # Try with _R1 suffix
        path_r1 = os.path.join(src_dir, sample_name + "_R1" + ext)
        if os.path.isfile(path_r1):
            candidates.append(path_r1)
        # Try with _1 suffix
        path_1 = os.path.join(src_dir, sample_name + "_1" + ext)
        if os.path.isfile(path_1):
            candidates.append(path_1)

    if candidates:
        return candidates

    # Glob fallback: any file starting with sample_name
    pattern = os.path.join(src_dir, sample_name + "*")
    matches = _glob.glob(pattern)
    # Filter to only fastq-like files and only R1/unpaired
    for m in sorted(matches):
        if os.path.isfile(m):
            lower = m.lower()
            if any(lower.endswith(e) for e in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]):
                # Skip R2/_2 to avoid double-linking
                base = os.path.basename(m)
                if "_R2" in base or "_2." in base:
                    continue
                candidates.append(m)

    return candidates


@app.get("/", response_class=HTMLResponse)
def root() -> str:
    return """
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <title>MONSDA Web Configurator</title>
  <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css" rel="stylesheet">
  <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap" rel="stylesheet">
  <style>
    * { box-sizing: border-box; }
    body { font-family: 'Inter', system-ui, sans-serif; margin: 0; padding: 32px; background: linear-gradient(135deg, #f0f4ff 0%, #fafbff 50%, #f5f3ff 100%); color: #1e293b; min-height: 100vh; }
    h1 { font-size: 1.8rem; font-weight: 700; margin-bottom: 4px; background: linear-gradient(135deg, #3b82f6, #8b5cf6); -webkit-background-clip: text; -webkit-text-fill-color: transparent; }
    h3 { margin: 0; font-size: 1.05rem; font-weight: 600; }
    .card { background: #fff; border: 1px solid #e2e8f0; border-radius: 16px; padding: 20px 24px; margin-bottom: 20px; box-shadow: 0 1px 3px rgba(0,0,0,.04), 0 4px 12px rgba(0,0,0,.02); transition: box-shadow .2s; }
    .card:hover { box-shadow: 0 4px 16px rgba(0,0,0,.06); }
    .row { display: grid; grid-template-columns: 1fr 1fr; gap: 16px; }
    label { font-size: .88rem; font-weight: 500; color: #475569; display: block; margin-bottom: 2px; }
    input, textarea, select { width: 100%; padding: 10px 14px; border: 1px solid #e2e8f0; border-radius: 10px; font-size: .9rem; margin-top: 4px; margin-bottom: 12px; transition: border-color .15s, box-shadow .15s; outline: none; }
    input:focus, textarea:focus { border-color: #3b82f6; box-shadow: 0 0 0 3px rgba(59,130,246,.1); }
    textarea { min-height: 140px; font-family: 'JetBrains Mono', ui-monospace, monospace; font-size: .82rem; }
    button { padding: 10px 18px; border: none; border-radius: 10px; font-size: .88rem; font-weight: 500; cursor: pointer; transition: all .15s; }
    .btn-primary, button[onclick*="preview"], button[onclick*="save"], button[onclick*="parse"], button[onclick*="create"] { background: linear-gradient(135deg, #3b82f6, #6366f1); color: #fff; box-shadow: 0 2px 8px rgba(59,130,246,.25); }
    .btn-primary:hover, button[onclick*="preview"]:hover, button[onclick*="save"]:hover, button[onclick*="parse"]:hover, button[onclick*="create"]:hover { transform: translateY(-1px); box-shadow: 0 4px 12px rgba(59,130,246,.35); }
    button[onclick*="loadTemplate"], button[onclick*="openToolsModal"], button[onclick*="openPathBrowser"] { background: #f1f5f9; color: #475569; border: 1px solid #e2e8f0; }
    button[onclick*="loadTemplate"]:hover, button[onclick*="openToolsModal"]:hover, button[onclick*="openPathBrowser"]:hover { background: #e2e8f0; }
    .ok { color: #059669; font-weight: 500; }
    .err { color: #dc2626; font-weight: 500; }
    .muted { color: #64748b; font-size: 0.85rem; }
    .pill { display: inline-block; padding: 3px 10px; border-radius: 999px; background: linear-gradient(135deg, #ede9fe, #e0e7ff); color: #4338ca; font-size: .82rem; font-weight: 500; margin-right: 6px; margin-bottom: 6px; }
    .pathline { display: flex; gap: 8px; align-items: center; margin-bottom: 12px; }
    .pathline input { flex: 1; margin-bottom: 0; }
    .pathline button { width: auto; white-space: nowrap; padding: 10px 14px; margin: 0; }
    .section-title { display: flex; align-items: center; gap: 10px; margin-bottom: 12px; }
    .help-btn { width: auto; margin: 0; padding: 3px 10px; border: 1px solid #e2e8f0; border-radius: 999px; background: #f8fafc; cursor: pointer; font-size: .78rem; color: #64748b; }
    .help-btn:hover { background: #e2e8f0; }
    .help-box { display: none; border: 1px solid #e2e8f0; background: linear-gradient(135deg, #f0f9ff, #faf5ff); border-radius: 10px; padding: 12px 16px; margin-bottom: 14px; font-size: .88rem; color: #475569; line-height: 1.6; }
    .quickstart { background: linear-gradient(135deg, #ecfeff, #f0f9ff); border-color: #bae6fd; }
    .workflow-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 8px; margin-top: 10px; }
    .wf-chip { display: flex; align-items: center; gap: 6px; border: 1px solid #e2e8f0; border-radius: 10px; padding: 8px 12px; background: #f8fafc; font-size: .88rem; cursor: pointer; transition: all .15s; }
    .wf-chip:hover { background: #ede9fe; border-color: #c4b5fd; }
    .wf-chip input { width: auto; margin: 0; }
    .inline-buttons { display: flex; gap: 10px; flex-wrap: wrap; margin-bottom: 10px; }
    .inline-buttons button { width: auto; }
    .cond-card { border: 1px solid #e2e8f0; border-radius: 12px; padding: 16px 20px; margin-bottom: 12px; background: #f8fafc; }
    .cond-card-header { display: flex; align-items: center; gap: 10px; margin-bottom: 12px; }
    .cond-card-header .cond-badge { background: linear-gradient(135deg, #6366f1, #8b5cf6); color: #fff; padding: 4px 12px; border-radius: 999px; font-size: .82rem; font-weight: 600; }
    .cond-card-header .cond-samples { color: #64748b; font-size: .82rem; }
    .cond-fields { display: grid; grid-template-columns: 1fr 1fr; gap: 10px 16px; }
    .cond-fields label { font-size: .82rem; margin-bottom: 0; }
    .cond-fields .pathline { margin-bottom: 0; }
    .cond-fields .pathline input { font-size: .82rem; padding: 7px 10px; margin-bottom: 0; }
    .cond-fields .pathline button { padding: 7px 10px; font-size: .78rem; }
    .cond-fields select { font-size: .82rem; padding: 7px 10px; margin-top: 4px; margin-bottom: 0; }
    .mode-toggle .mode-opt { background: #f8fafc; color: #64748b; }
    .mode-toggle .mode-opt.active { background: linear-gradient(135deg, #3b82f6, #6366f1); color: #fff; font-weight: 600; }
    .cond-file-list { list-style: none; padding: 0; margin: 4px 0 0 0; }
    .cond-file-list li { display: flex; align-items: center; gap: 6px; padding: 3px 0; font-size: .82rem; color: #334155; }
    .cond-file-list li button { background: none; border: none; color: #dc2626; cursor: pointer; padding: 0 4px; font-size: 1rem; }
    .overlay-backdrop { position: fixed; inset: 0; background: rgba(17, 24, 39, 0.55); display: none; align-items: center; justify-content: center; z-index: 1000; }
    .overlay-modal { background: #fff; width: min(900px, 94vw); max-height: 84vh; overflow: auto; border-radius: 10px; border: 1px solid #d1d5db; padding: 14px; }
    .fb-nav-btn { background:#fff; border:1px solid #e2e8f0; border-radius:8px; width:32px; height:32px; display:flex; align-items:center; justify-content:center; cursor:pointer; transition:all .15s; }
    .fb-nav-btn:hover { background:#e2e8f0; }
    .fb-crumb { color:#64748b; text-decoration:none; padding:2px 6px; border-radius:4px; transition:background .12s; white-space:nowrap; }
    .fb-crumb:hover { background:#e2e8f0; color:#1e293b; text-decoration:none; }
    .fb-crumb-sep { color:#cbd5e1; margin:0 1px; }
    .fb-crumb-active { color:#1e293b; font-weight:600; padding:2px 6px; }
    .fb-root-btn { display:block; padding:7px 16px; border:none; background:none; width:100%; text-align:left; font-size:.88rem; color:#475569; cursor:pointer; border-radius:0; transition:background .12s; }
    .fb-root-btn:hover { background:#e2e8f0; }
    .fb-entry { display:flex; align-items:center; gap:10px; padding:8px 20px; cursor:pointer; border:none; background:none; width:100%; text-align:left; font-size:.9rem; color:#334155; transition:background .1s; user-select:none; }
    .fb-entry:hover { background:#f1f5f9; }
    .fb-entry.selected { background:#dbeafe; }
    .fb-entry-icon { flex-shrink:0; width:20px; height:20px; display:flex; align-items:center; justify-content:center; }
    .fb-entry-icon.folder { color:#f59e0b; }
    .fb-entry-icon.file { color:#64748b; }
    .fb-entry-name { flex:1; overflow:hidden; text-overflow:ellipsis; white-space:nowrap; }
    .tools-group { border: 1px solid #e2e8f0; border-radius: 12px; padding: 12px 16px; margin-bottom: 10px; background: #f8fafc; }
    .tools-group summary { cursor: pointer; font-weight: 600; color: #334155; }
    .tools-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 6px 12px; margin-top: 10px; }
    .tools-grid label { font-size: 0.88rem; display: flex; align-items: center; gap: 6px; padding: 4px 8px; border-radius: 6px; transition: background .15s; }
    .tools-grid label:hover { background: #ede9fe; }
  </style>
</head>
<body>
  <h1>MONSDA Web Configurator</h1>
  <p class="muted">Interactive builder for MONSDA pipeline configurations.</p>

  <div class="card quickstart">
    <div class="section-title">
      <h3>Quick Start</h3>
      <button type="button" class="help-btn" onclick="toggleHelp('helpQuick')">?</button>
    </div>
    <div id="helpQuick" class="help-box" style="display:block;">
      <strong>1.</strong> Pick and parse a samplesheet.<br/>
      <strong>2.</strong> Select workflows &amp; configure tools.<br/>
      <strong>3.</strong> Preview and save config.<br/>
      <strong>4.</strong> Create project skeleton with symlinked data.
    </div>
  </div>

  <div class="card">
    <div class="section-title">
      <h3>1. Samplesheet &rarr; Settings</h3>
      <button type="button" class="help-btn" onclick="toggleHelp('helpSamplesheet')">?</button>
    </div>
    <div id="helpSamplesheet" class="help-box">
      Choose your CSV/TSV samplesheet, then click <b>Parse Samplesheet</b>. This fills the SETTINGS JSON automatically.
    </div>
    <label>Samplesheet path (CSV/TSV)</label>
    <div class="pathline">
      <input id="samplesheetPath" placeholder="/abs/path/to/samplesheet.csv" />
      <button type="button" onclick="openPathBrowser('samplesheetPath','samplesheets')">Browse…</button>
    </div>
    <button type="button" onclick="parseSamplesheet()">Parse Samplesheet</button>
    <div id="samplesheetStatus"></div>
    <div id="conditionsPreview"></div>
  </div>

  <div class="card">
    <div class="section-title">
      <h3>2. Workflows &amp; Tools</h3>
      <button type="button" class="help-btn" onclick="toggleHelp('helpWorkflows')">?</button>
    </div>
    <div id="helpWorkflows" class="help-box">
      Select the workflows you need. Tools are only configured for workflows you select. Click <b>Configure tools</b> to customize per workflow.
    </div>
    <div class="inline-buttons">
      <button type="button" onclick="loadTemplate()">Reload workflows from template</button>
      <button type="button" onclick="openToolsModal()">Configure tools for selected workflows</button>
    </div>
    <div id="templateStatus" class="muted"></div>
    <div id="workflowSummary" class="muted"></div>
    <div id="workflowChooser" class="workflow-grid"></div>
  </div>

  <div class="card">
    <div class="section-title">
      <h3>3. Output &amp; Finalize</h3>
      <button type="button" class="help-btn" onclick="toggleHelp('helpBuild')">?</button>
    </div>
    <div id="helpBuild" class="help-box">
      Choose one mode:<br/>
      <strong>Save config only</strong> &mdash; generates and saves the JSON config to a directory of your choice.<br/>
      <strong>Create project</strong> &mdash; builds a full MONSDA project skeleton with FASTQ &amp; GENOMES symlinks, plus the config. Select FASTQ files per condition either by directory (filtered by sample names) or by picking individual files.
    </div>

    <!-- Mode toggle -->
    <div class="mode-toggle" style="display:flex; gap:0; margin-bottom:16px; border:1px solid #e2e8f0; border-radius:10px; overflow:hidden;">
      <label class="mode-opt" style="flex:1; text-align:center; padding:10px; cursor:pointer; font-size:.88rem; font-weight:500; transition:all .15s;" id="modeOptConfig">
        <input type="radio" name="outputMode" value="config" checked onchange="switchOutputMode('config')" style="display:none;"/>
        Save config only
      </label>
      <label class="mode-opt" style="flex:1; text-align:center; padding:10px; cursor:pointer; font-size:.88rem; font-weight:500; transition:all .15s;" id="modeOptProject">
        <input type="radio" name="outputMode" value="project" onchange="switchOutputMode('project')" style="display:none;"/>
        Create project
      </label>
    </div>

    <!-- Shared fields -->
    <div class="row">
      <div>
        <label>Config / project name</label>
        <input id="configName" value="monsda" />
      </div>
      <div>
        <label>MAXTHREADS</label>
        <input id="maxthreads" value="16" />
      </div>
    </div>

    <label>SETTINGS JSON (auto-filled from samplesheet, editable)</label>
    <textarea id="settingsJson"></textarea>

    <!-- Config-only panel -->
    <div id="panelConfigOnly">
      <label>Output directory</label>
      <div class="pathline">
        <input id="outputDir" placeholder="/abs/output/dir" />
        <button type="button" onclick="openPathBrowser('outputDir','dirs')">Browse&hellip;</button>
      </div>
      <div class="inline-buttons">
        <button type="button" onclick="previewConfig()">Preview Config JSON</button>
        <button type="button" onclick="saveConfig()">Save Config</button>
      </div>
    </div>

    <!-- Project panel -->
    <div id="panelProject" style="display:none;">
      <div class="row">
        <div>
          <label>Project directory (parent folder)</label>
          <div class="pathline">
            <input id="projectDir" placeholder="/abs/path/to/parent" />
            <button type="button" onclick="openPathBrowser('projectDir','dirs')">Browse&hellip;</button>
          </div>
        </div>
        <div>
          <label class="muted" style="margin-top:22px;">A subfolder with the config name above will be created here.</label>
        </div>
      </div>

      <div class="inline-buttons" style="margin-top:8px;">
        <button type="button" onclick="loadProjectConditions()">Load conditions from settings</button>
      </div>

      <div id="conditionWizard" style="margin-top:16px;"></div>

      <div class="inline-buttons" style="margin-top:16px;">
        <button type="button" onclick="previewConfig()">Preview Config JSON</button>
        <button type="button" onclick="createProject()">Create Project</button>
      </div>
    </div>

    <div id="buildStatus"></div>
    <div id="projectStatus"></div>

    <label>Config preview</label>
    <textarea id="configPreview"></textarea>
  </div>

  <div class="modal fade" id="browserModal" tabindex="-1" aria-hidden="true">
    <div class="modal-dialog modal-xl modal-dialog-centered" style="max-width:860px;">
      <div class="modal-content" style="border-radius:16px; overflow:hidden; box-shadow: 0 25px 60px rgba(0,0,0,.25);">
        <!-- Header -->
        <div style="background:linear-gradient(135deg,#1e293b,#334155); color:#fff; padding:16px 24px; display:flex; align-items:center; gap:12px;">
          <svg width="22" height="22" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M22 19a2 2 0 0 1-2 2H4a2 2 0 0 1-2-2V5a2 2 0 0 1 2-2h5l2 3h9a2 2 0 0 1 2 2z"/></svg>
          <span style="font-size:1.1rem; font-weight:600;" id="browserTitle">Select a path</span>
          <span style="flex:1;"></span>
          <button type="button" style="background:none; border:none; color:#94a3b8; font-size:1.4rem; cursor:pointer;" data-bs-dismiss="modal" aria-label="Close">&times;</button>
        </div>
        <!-- Location bar -->
        <div style="background:#f8fafc; border-bottom:1px solid #e2e8f0; padding:10px 20px; display:flex; gap:8px; align-items:center;">
          <button type="button" class="fb-nav-btn" onclick="browseUp()" title="Go up">
            <svg width="18" height="18" fill="none" stroke="currentColor" stroke-width="2"><path d="M15 11l-6-6-6 6"/></svg>
          </button>
          <div id="browseBreadcrumb" style="display:flex; align-items:center; gap:2px; flex:1; overflow-x:auto; font-size:.9rem;"></div>
        </div>
        <!-- Body: sidebar + file list -->
        <div style="display:flex; min-height:420px; max-height:65vh;">
          <!-- Sidebar -->
          <div style="width:180px; background:#f1f5f9; border-right:1px solid #e2e8f0; padding:12px 0; overflow-y:auto; flex-shrink:0;">
            <div style="padding:0 12px; margin-bottom:8px; font-size:.75rem; text-transform:uppercase; color:#64748b; font-weight:600; letter-spacing:.03em;">Locations</div>
            <div id="browseRoots"></div>
          </div>
          <!-- Main file list -->
          <div style="flex:1; overflow-y:auto; padding:8px 0;" id="browseListContainer">
            <div id="browseList"></div>
            <div id="browseEmpty" style="display:none; padding:48px 24px; text-align:center; color:#94a3b8;">
              <svg width="48" height="48" fill="none" stroke="currentColor" stroke-width="1.5" style="margin:0 auto 12px;"><path d="M3 7v30a4 4 0 0 0 4 4h34a4 4 0 0 0 4-4V15a4 4 0 0 0-4-4H24l-4-4H7a4 4 0 0 0-4 4z"/></svg>
              <div>This folder is empty</div>
            </div>
          </div>
        </div>
        <!-- Footer -->
        <div style="background:#f8fafc; border-top:1px solid #e2e8f0; padding:12px 20px; display:flex; align-items:center; gap:10px;">
          <input id="browsePathInput" class="form-control form-control-sm" style="flex:1; border-radius:8px;" placeholder="Path..." />
          <button type="button" class="btn btn-sm btn-outline-secondary" style="border-radius:8px;" onclick="browseTo(document.getElementById('browsePathInput').value)">Go</button>
          <button type="button" class="btn btn-sm btn-primary" style="border-radius:8px; min-width:120px;" onclick="chooseCurrentPath()">Select</button>
        </div>
      </div>
    </div>
  </div>

  <div id="toolsModal" class="overlay-backdrop" onclick="if(event.target===this) closeToolsModal()">
    <div class="overlay-modal">
      <div class="section-title">
        <h3>Tool selection for active workflows</h3>
        <button type="button" class="help-btn" onclick="closeToolsModal()">Close</button>
      </div>
      <div class="muted">By default all tools are selected when you enable a workflow. Uncheck tools you do not want.</div>
      <div id="toolsModalBody"></div>
    </div>
  </div>

<script>
let templateFields = null;
let lastConfig = null;
let browserState = { target: '', mode: 'dirs', current: '' };
let selectedToolsByWorkflow = {};
let browserModalRef = null;
let selectedBrowsePath = '';
let currentOutputMode = 'config';

function switchOutputMode(mode) {
  currentOutputMode = mode;
  document.getElementById('panelConfigOnly').style.display = mode === 'config' ? '' : 'none';
  document.getElementById('panelProject').style.display = mode === 'project' ? '' : 'none';
  document.getElementById('modeOptConfig').classList.toggle('active', mode === 'config');
  document.getElementById('modeOptProject').classList.toggle('active', mode === 'project');
}

function toggleHelp(id) {
  const el = document.getElementById(id);
  if (!el) return;
  el.style.display = (el.style.display === 'block') ? 'none' : 'block';
}

function setStatus(id, text, ok=true) {
  const el = document.getElementById(id);
  if (!el) return;
  el.textContent = text;
  el.className = ok ? 'ok' : 'err';
}

function escHtml(s) {
  return String(s)
    .replaceAll('&', '&amp;')
    .replaceAll('<', '&lt;')
    .replaceAll('>', '&gt;')
    .replaceAll('"', '&quot;')
    .replaceAll("'", '&#039;');
}

function getToolsForWorkflow(wf) {
  if (!templateFields || !templateFields[wf] || !templateFields[wf].TOOLS) return [];
  return Object.keys(templateFields[wf].TOOLS);
}

function onWorkflowToggle(wf, checked) {
  if (checked && !selectedToolsByWorkflow[wf]) {
    selectedToolsByWorkflow[wf] = [];
  }
  const selected = Array.from(document.querySelectorAll('input[data-workflow]:checked')).map(i => i.getAttribute('data-workflow'));
  document.getElementById('workflowSummary').textContent = selected.length ? `Selected workflows: ${selected.join(', ')}` : 'No workflows selected yet.';
}

async function loadTemplate() {
  const holder = document.getElementById('workflowChooser');
  const status = document.getElementById('templateStatus');
  const summary = document.getElementById('workflowSummary');
  try {
    status.textContent = 'Loading template...';
    const r = await fetch('/template/fields');
    const data = await r.json();
    if (!r.ok) {
      throw new Error(data.detail || 'Failed to load template fields');
    }

    templateFields = data;
    holder.innerHTML = '';

    const workflows = Object.keys(data)
      .filter(k => !['WORKFLOWS','BINS','MAXTHREADS','SETTINGS','VERSION'].includes(k))
      .sort();

    workflows.forEach(wf => {
      holder.innerHTML += `<label class="wf-chip"><input type="checkbox" data-workflow="${wf}" onchange="onWorkflowToggle('${wf}', this.checked)" />${escHtml(wf)}</label>`;
    });

    status.textContent = `Loaded ${workflows.length} workflows from template.`;
    status.className = 'ok';
    summary.textContent = 'No workflows selected yet.';
  } catch (e) {
    holder.innerHTML = '';
    summary.textContent = '';
    status.textContent = `Failed to load template fields: ${e.message || e}`;
    status.className = 'err';
  }
}

function openToolsModal() {
  const modal = document.getElementById('toolsModal');
  const body = document.getElementById('toolsModalBody');
  const selected = Array.from(document.querySelectorAll('input[data-workflow]:checked')).map(i => i.getAttribute('data-workflow'));
  if (!selected.length) {
    body.innerHTML = '<div class="muted">Select at least one workflow first.</div>';
    modal.style.display = 'flex';
    return;
  }

  body.innerHTML = selected.map(wf => {
    const tools = getToolsForWorkflow(wf);
    const active = new Set(selectedToolsByWorkflow[wf] || []);
    const toolHtml = tools.length
      ? tools.map(t => `<label><input type="checkbox" ${active.has(t) ? 'checked' : ''} onchange="onToolToggle('${wf}','${t}',this.checked)"/> ${escHtml(t)}</label>`).join('')
      : '<div class="muted">No tool list available.</div>';
    return `<details class="tools-group" open><summary>${escHtml(wf)}</summary><div class="tools-grid">${toolHtml}</div></details>`;
  }).join('');
  modal.style.display = 'flex';
}

function closeToolsModal() {
  document.getElementById('toolsModal').style.display = 'none';
}

function onToolToggle(wf, tool, checked) {
  const current = new Set(selectedToolsByWorkflow[wf] || []);
  if (checked) current.add(tool);
  else current.delete(tool);
  selectedToolsByWorkflow[wf] = Array.from(current);
}

async function loadBrowseRoots() {
  const r = await fetch('/fs/roots');
  const data = await r.json();
  if (!r.ok) throw new Error(data.detail || 'Failed to load roots');
  const el = document.getElementById('browseRoots');
  el.innerHTML = '';
  (data.roots || []).forEach(p => {
    const label = p === '/' ? 'Root (/)' : p.split('/').filter(Boolean).pop() || p;
    const btn = document.createElement('button');
    btn.type = 'button';
    btn.className = 'fb-root-btn';
    btn.innerHTML = `<svg width="14" height="14" fill="none" stroke="currentColor" stroke-width="2" style="margin-right:6px; vertical-align:-2px;"><path d="M3 3h4l1.5 2H13a1 1 0 0 1 1 1v6a1 1 0 0 1-1 1H3a1 1 0 0 1-1-1V4a1 1 0 0 1 1-1z"/></svg>${escHtml(label)}`;
    btn.addEventListener('click', () => browseTo(p));
    el.appendChild(btn);
  });
}

function renderBreadcrumb(path) {
  const crumb = document.getElementById('browseBreadcrumb');
  if (!crumb) return;
  crumb.innerHTML = '';
  const parts = (path || '/').split('/').filter(Boolean);
  let acc = '/';

  const rootLink = document.createElement('a');
  rootLink.href = '#';
  rootLink.className = 'fb-crumb';
  rootLink.textContent = '/';
  rootLink.addEventListener('click', (e) => { e.preventDefault(); browseTo('/'); });
  crumb.appendChild(rootLink);

  for (let i = 0; i < parts.length; i++) {
    const p = parts[i];
    acc = acc === '/' ? `/${p}` : `${acc}/${p}`;

    const sep = document.createElement('span');
    sep.className = 'fb-crumb-sep';
    sep.textContent = '/';
    crumb.appendChild(sep);

    if (i === parts.length - 1) {
      const span = document.createElement('span');
      span.className = 'fb-crumb-active';
      span.textContent = p;
      crumb.appendChild(span);
    } else {
      const link = document.createElement('a');
      link.href = '#';
      link.className = 'fb-crumb';
      link.textContent = p;
      const target = acc;
      link.addEventListener('click', (e) => { e.preventDefault(); browseTo(target); });
      crumb.appendChild(link);
    }
  }
}

async function browseTo(path) {
  try {
    const q = new URLSearchParams();
    q.set('mode', browserState.mode || 'dirs');
    if (path) q.set('path', path);
    const r = await fetch(`/fs/list?${q.toString()}`);
    const data = await r.json();
    if (!r.ok) throw new Error(data.detail || 'Browse failed');

    browserState.current = data.current;
    selectedBrowsePath = data.current;
    document.getElementById('browsePathInput').value = data.current;
    renderBreadcrumb(data.current);

    const modeLabel = browserState.mode === 'dirs' ? 'folder' : 'file';
    document.getElementById('browserTitle').textContent = `Select a ${modeLabel}`;

    const folderSvg = `<svg width="20" height="20" fill="none" stroke="currentColor" stroke-width="1.5"><path d="M2 5v11a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V9a2 2 0 0 0-2-2h-6l-2-2H4a2 2 0 0 0-2 2z"/></svg>`;
    const fileSvg = `<svg width="20" height="20" fill="none" stroke="currentColor" stroke-width="1.5"><path d="M14 2H6a2 2 0 0 0-2 2v12a2 2 0 0 0 2 2h8a2 2 0 0 0 2-2V6l-4-4z"/><polyline points="14 2 14 6 18 6"/></svg>`;

    const list = document.getElementById('browseList');
    const empty = document.getElementById('browseEmpty');
    list.innerHTML = '';

    (data.dirs || []).forEach(d => {
      const btn = document.createElement('button');
      btn.type = 'button';
      btn.className = 'fb-entry';
      btn.innerHTML = `<span class="fb-entry-icon folder">${folderSvg}</span><span class="fb-entry-name">${escHtml(d.name)}</span>`;
      btn.addEventListener('click', () => markBrowseSelection(d.path, btn));
      btn.addEventListener('dblclick', () => browseTo(d.path));
      list.appendChild(btn);
    });

    (data.files || []).forEach(f => {
      const btn = document.createElement('button');
      btn.type = 'button';
      btn.className = 'fb-entry';
      btn.innerHTML = `<span class="fb-entry-icon file">${fileSvg}</span><span class="fb-entry-name">${escHtml(f.name)}</span>`;
      btn.addEventListener('click', () => { markBrowseSelection(f.path, btn); choosePath(f.path); });
      list.appendChild(btn);
    });

    empty.style.display = list.children.length ? 'none' : 'block';
  } catch (e) {
    setStatus('browseStatus', e.message || String(e), false);
  }
}

function markBrowseSelection(path, el) {
  selectedBrowsePath = path;
  document.getElementById('browsePathInput').value = path;
  const list = document.getElementById('browseList');
  if (!list) return;
  list.querySelectorAll('.fb-entry').forEach(x => x.classList.remove('selected'));
  if (el) el.classList.add('selected');
}

function openPathBrowser(target, mode='dirs') {
  browserState.target = target;
  browserState.mode = mode;
  if (!browserModalRef) {
    const el = document.getElementById('browserModal');
    browserModalRef = new bootstrap.Modal(el, { backdrop: true, keyboard: true });
  }
  browserModalRef.show();
  const start = (document.getElementById(target)?.value || '').trim();
  loadBrowseRoots().then(() => browseTo(start || '')).catch((e) => {
    setStatus('browseStatus', e.message || String(e), false);
  });
}

function closePathBrowser() {
  if (browserModalRef) browserModalRef.hide();
}

function choosePath(path) {
  if (!browserState.target) return;
  document.getElementById(browserState.target).value = path;
  setStatus('browseStatus', `Selected: ${path}`, true);
  closePathBrowser();
}

function chooseCurrentPath() {
  const picked = selectedBrowsePath || browserState.current;
  if (!picked) return;
  choosePath(picked);
}

function browseUp() {
  if (!browserState.current) return;
  const cur = browserState.current;
  if (cur === '/') {
    browseTo('/');
    return;
  }
  const parts = cur.split('/').filter(Boolean);
  parts.pop();
  const parent = '/' + parts.join('/');
  browseTo(parent || '/');
}

async function parseSamplesheet() {
  const p = document.getElementById('samplesheetPath').value.trim();
  if (!p) {
    setStatus('samplesheetStatus', 'Please provide a samplesheet path.', false);
    return;
  }
  try {
    const r = await fetch('/samplesheet/parse', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ samplesheet_path: p })
    });
    const data = await r.json();
    if (r.ok) {
      document.getElementById('settingsJson').value = JSON.stringify(data.settings, null, 2);
      document.getElementById('conditionsPreview').innerHTML = data.conditions.map(c => `<span class="pill">${c}</span>`).join('');
      setStatus('samplesheetStatus', `Parsed ${data.conditions.length} condition(s) from samplesheet.`);
    } else {
      setStatus('samplesheetStatus', data.detail || 'Failed to parse samplesheet.', false);
    }
  } catch (e) {
    setStatus('samplesheetStatus', 'Failed to parse samplesheet.', false);
  }
}

function collectWorkflowSelection() {
  const selectedWorkflows = Array.from(document.querySelectorAll('input[data-workflow]:checked')).map(i => i.getAttribute('data-workflow'));
  const tools = {};
  selectedWorkflows.forEach(wf => {
    const allTools = getToolsForWorkflow(wf);
    const picked = selectedToolsByWorkflow[wf] || allTools;
    tools[wf] = picked.length ? picked : allTools;
  });
  return { selectedWorkflows, tools };
}

function parseSettingsJson() {
  const txt = document.getElementById('settingsJson').value.trim();
  if (!txt) return null;
  return JSON.parse(txt);
}

async function previewConfig() {
  try {
    const { selectedWorkflows, tools } = collectWorkflowSelection();
    const settings = parseSettingsJson();
    const body = {
      config_name: document.getElementById('configName').value.trim(),
      output_dir: document.getElementById('outputDir').value.trim(),
      workflows: selectedWorkflows,
      tools: tools,
      maxthreads: document.getElementById('maxthreads').value.trim(),
      settings: settings,
      samplesheet_path: document.getElementById('samplesheetPath').value.trim() || null
    };

    const r = await fetch('/config/preview', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(body)
    });
    const data = await r.json();
    if (r.ok) {
      lastConfig = data.config;
      document.getElementById('configPreview').value = JSON.stringify(data.config, null, 2);
      setStatus('buildStatus', 'Config preview generated.');
    } else {
      setStatus('buildStatus', data.detail || 'Failed to preview config.', false);
    }
  } catch (e) {
    setStatus('buildStatus', 'Invalid SETTINGS JSON or request failed.', false);
  }
}

async function saveConfig() {
  try {
    if (!lastConfig) {
      setStatus('buildStatus', 'Generate preview first.', false);
      return;
    }
    const r = await fetch('/config/save', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        config_name: document.getElementById('configName').value.trim(),
        output_dir: document.getElementById('outputDir').value.trim(),
        config: lastConfig
      })
    });
    const data = await r.json();
    if (r.ok) {
      setStatus('buildStatus', `Saved: ${data.path}`);
    } else {
      setStatus('buildStatus', data.detail || 'Failed to save config.', false);
    }
  } catch (e) {
    setStatus('buildStatus', 'Failed to save config.', false);
  }
}

async function createProject() {
  try {
    const projectDir = document.getElementById('projectDir').value.trim();
    const projectName = document.getElementById('configName').value.trim() || 'monsda';
    if (!projectDir) {
      setStatus('projectStatus', 'Please set a project directory.', false);
      return;
    }

    // Collect per-condition file selections from wizard
    const conditionFiles = [];
    document.querySelectorAll('.cond-card').forEach(card => {
      const cond = card.dataset.condition;
      // Determine FASTQ mode
      const modeRadio = card.querySelector('input[type="radio"][value="files"]:checked');
      const useFiles = !!modeRadio;
      let fastqFiles = [];
      if (useFiles) {
        card.querySelectorAll('.cond-file-list li').forEach(li => {
          if (li.dataset.path) fastqFiles.push(li.dataset.path);
        });
      }
      conditionFiles.push({
        condition: cond,
        fastq_dir: useFiles ? '' : (card.querySelector('.cond-fastq-dir')?.value?.trim() || ''),
        fastq_files: fastqFiles,
        sequencing: card.querySelector('.cond-sequencing')?.value || 'paired',
        reference: card.querySelector('.cond-reference')?.value?.trim() || '',
        gtf: card.querySelector('.cond-gtf')?.value?.trim() || '',
        gff: card.querySelector('.cond-gff')?.value?.trim() || '',
        decoy: card.querySelector('.cond-decoy')?.value?.trim() || ''
      });
    });

    const r = await fetch('/project/create', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        project_dir: projectDir,
        project_name: projectName,
        condition_files: conditionFiles,
        config: lastConfig || null,
        settings: parseSettingsJson()
      })
    });
    const data = await r.json();
    if (r.ok) {
      let msg = `Project created at <strong>${escHtml(data.path)}</strong>. ${data.linked_files} file(s) symlinked.`;
      if (data.config_path) msg += `<br/>Config written: ${escHtml(data.config_path)}`;
      if (data.warnings && data.warnings.length) msg += `<br/><span class="err">Warnings:</span> ${data.warnings.map(escHtml).join('<br/>')}`;
      document.getElementById('projectStatus').innerHTML = `<span class="ok">${msg}</span>`;
    } else {
      setStatus('projectStatus', data.detail || 'Failed to create project.', false);
    }
  } catch (e) {
    setStatus('projectStatus', 'Failed to create project: ' + e.message, false);
  }
}

function loadProjectConditions() {
  const settings = parseSettingsJson();
  if (!settings) {
    setStatus('projectStatus', 'Parse a samplesheet first or edit SETTINGS JSON manually, then load conditions.', false);
    return;
  }

  // Extract condition leaves from settings (those with SAMPLES key)
  const conditions = [];
  function walkSettings(node, path) {
    if (!node || typeof node !== 'object') return;
    if (Array.isArray(node.SAMPLES)) {
      conditions.push({ path: path.join('/'), samples: node.SAMPLES, sequencing: node.SEQUENCING || 'paired', reference: node.REFERENCE || '', gtf: (node.ANNOTATION && node.ANNOTATION.GTF) || '', gff: (node.ANNOTATION && node.ANNOTATION.GFF) || '', decoy: node.DECOY || '' });
      return;
    }
    for (const [k, v] of Object.entries(node)) {
      if (typeof v === 'object' && !Array.isArray(v)) walkSettings(v, [...path, k]);
    }
  }
  walkSettings(settings, []);

  if (!conditions.length) {
    setStatus('projectStatus', 'No conditions with SAMPLES found in settings.', false);
    return;
  }

  const wizard = document.getElementById('conditionWizard');
  wizard.innerHTML = `<div class="muted" style="margin-bottom:10px;">${conditions.length} condition(s) found. For each, select FASTQ files by directory (auto-matched by sample name) or pick individual files.</div>`;

  conditions.forEach((c, idx) => {
    const id = `cond_${idx}`;
    const hasSamples = c.samples && c.samples.length > 0;
    const samplesStr = hasSamples
      ? c.samples.slice(0, 5).join(', ') + (c.samples.length > 5 ? ` (+${c.samples.length - 5} more)` : '')
      : '(no samples defined)';
    const card = document.createElement('div');
    card.className = 'cond-card';
    card.dataset.condition = c.path;
    card.innerHTML = `
      <div class="cond-card-header">
        <span class="cond-badge">${escHtml(c.path)}</span>
        <span class="cond-samples">${escHtml(samplesStr)}</span>
      </div>
      <div class="cond-fields">
        <div style="grid-column:1/-1;">
          <label>FASTQ source mode</label>
          <div style="display:flex; gap:12px; align-items:center; margin-top:4px;">
            <label style="font-weight:400; display:inline-flex; align-items:center; gap:4px; margin:0;">
              <input type="radio" name="${id}_fqmode" value="dir" checked onchange="toggleFqMode('${id}','dir')" style="width:auto; margin:0;"/> Directory (filter by sample names)
            </label>
            <label style="font-weight:400; display:inline-flex; align-items:center; gap:4px; margin:0;">
              <input type="radio" name="${id}_fqmode" value="files" onchange="toggleFqMode('${id}','files')" style="width:auto; margin:0;"/> Pick individual files
            </label>
          </div>
        </div>
        <div class="fq-dir-panel" id="${id}_dirpanel">
          <label>FASTQ source directory</label>
          <div class="pathline">
            <input class="cond-fastq-dir" id="${id}_fastq" placeholder="/path/to/fastq/files" />
            <button type="button" onclick="openPathBrowser('${id}_fastq','dirs')">Browse&hellip;</button>
          </div>
        </div>
        <div class="fq-files-panel" id="${id}_filespanel" style="display:none;">
          <label>Selected FASTQ files <button type="button" style="width:auto; padding:2px 8px; font-size:.78rem; margin:0 0 0 8px;" onclick="addFastqFile('${id}')">+ Add file</button></label>
          <ul class="cond-file-list" id="${id}_filelist"></ul>
        </div>
        <div>
          <label>Sequencing</label>
          <select class="cond-sequencing">
            <option value="paired" ${c.sequencing === 'paired' ? 'selected' : ''}>Paired-end</option>
            <option value="single" ${c.sequencing !== 'paired' ? 'selected' : ''}>Single-end</option>
          </select>
        </div>
        <div>
          <label>Reference genome (.fa.gz)</label>
          <div class="pathline">
            <input class="cond-reference" id="${id}_ref" value="${escHtml(c.reference)}" placeholder="/path/to/genome.fa.gz" />
            <button type="button" onclick="openPathBrowser('${id}_ref','all')">Browse&hellip;</button>
          </div>
        </div>
        <div>
          <label>GTF annotation (.gtf.gz)</label>
          <div class="pathline">
            <input class="cond-gtf" id="${id}_gtf" value="${escHtml(c.gtf)}" placeholder="/path/to/annotation.gtf.gz" />
            <button type="button" onclick="openPathBrowser('${id}_gtf','all')">Browse&hellip;</button>
          </div>
        </div>
        <div>
          <label>GFF annotation (optional)</label>
          <div class="pathline">
            <input class="cond-gff" id="${id}_gff" value="${escHtml(c.gff)}" placeholder="(optional)" />
            <button type="button" onclick="openPathBrowser('${id}_gff','all')">Browse&hellip;</button>
          </div>
        </div>
        <div>
          <label>Decoy (optional)</label>
          <div class="pathline">
            <input class="cond-decoy" id="${id}_decoy" value="${escHtml(c.decoy)}" placeholder="(optional)" />
            <button type="button" onclick="openPathBrowser('${id}_decoy','all')">Browse&hellip;</button>
          </div>
        </div>
      </div>
    `;
    wizard.appendChild(card);
  });
  setStatus('projectStatus', `Loaded ${conditions.length} condition(s). Fill in file paths and click Create Project.`);
}

function toggleFqMode(id, mode) {
  document.getElementById(id + '_dirpanel').style.display = mode === 'dir' ? '' : 'none';
  document.getElementById(id + '_filespanel').style.display = mode === 'files' ? '' : 'none';
}

function addFastqFile(id) {
  // Use the file browser targeting a hidden input, then on selection add to the list
  const tempId = id + '_fqtemp';
  let tempInput = document.getElementById(tempId);
  if (!tempInput) {
    tempInput = document.createElement('input');
    tempInput.type = 'hidden';
    tempInput.id = tempId;
    document.body.appendChild(tempInput);
  }
  tempInput.value = '';
  // Override choosePath temporarily to add to file list instead
  const origChoose = choosePath;
  choosePath = function(path) {
    if (!browserState.target) return;
    // Add to file list
    const list = document.getElementById(id + '_filelist');
    const li = document.createElement('li');
    li.innerHTML = `<span>${escHtml(path)}</span>`;
    const rmBtn = document.createElement('button');
    rmBtn.type = 'button';
    rmBtn.textContent = '\u00d7';
    rmBtn.addEventListener('click', () => li.remove());
    li.appendChild(rmBtn);
    li.dataset.path = path;
    list.appendChild(li);
    closePathBrowser();
    choosePath = origChoose;
  };
  openPathBrowser(tempId, 'all');
}

window.onload = () => {
  loadTemplate();
  switchOutputMode('config');
};
</script>
<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/js/bootstrap.bundle.min.js"></script>
</body>
</html>
"""
