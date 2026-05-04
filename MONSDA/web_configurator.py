import copy
import json
import os
from typing import Any, Dict, List, Optional

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse
from pydantic import BaseModel, Field
from snakemake.common.configfile import load_configfile

from .Params import samplesheet_to_settings

TEMPLATE_PATH = os.path.join(
    os.path.dirname(__file__), "../configs/template_base_commented.json"
)
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


class ProjectRequest(BaseModel):
    project_dir: str
    workflows: List[str] = Field(default_factory=list)


class SaveConfigRequest(BaseModel):
    config_name: str
    output_dir: str
    config: Dict[str, Any]


def load_template() -> Dict[str, Any]:
    return load_configfile(TEMPLATE_PATH)


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
        if e.is_dir(follow_symlinks=False):
            dirs.append({"name": e.name, "path": e.path})
            continue
        if not e.is_file(follow_symlinks=False):
            continue
        if mode in {"all", "samplesheets"}:
            if mode == "samplesheets":
                lower = e.name.lower()
                if not lower.endswith((".csv", ".tsv", ".txt")):
                    continue
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
    project_dir = os.path.abspath(req.project_dir.strip())
    os.makedirs(project_dir, exist_ok=True)

    # Basic MONSDA project skeleton
    os.makedirs(os.path.join(project_dir, "FASTQ"), exist_ok=True)
    os.makedirs(os.path.join(project_dir, "GENOMES"), exist_ok=True)
    os.makedirs(os.path.join(project_dir, "LOGS"), exist_ok=True)

    for wf in req.workflows:
        os.makedirs(os.path.join(project_dir, wf), exist_ok=True)

    return {"message": "Project directory prepared.", "path": project_dir}


@app.get("/", response_class=HTMLResponse)
def root() -> str:
    return """
<!DOCTYPE html>
<html lang=\"en\">
<head>
  <meta charset=\"UTF-8\" />
  <title>MONSDA Web Configurator</title>
  <style>
    body { font-family: system-ui, sans-serif; margin: 24px; background: #f7f8fb; color: #1f2937; }
    h1 { margin-bottom: 8px; }
    .card { background: #fff; border: 1px solid #e5e7eb; border-radius: 10px; padding: 16px; margin-bottom: 16px; }
    .row { display: grid; grid-template-columns: 1fr 1fr; gap: 12px; }
    input, textarea, select, button { width: 100%; box-sizing: border-box; margin-top: 6px; margin-bottom: 10px; }
    textarea { min-height: 140px; font-family: ui-monospace, monospace; }
    .ok { color: #065f46; }
    .err { color: #991b1b; }
    .muted { color: #6b7280; font-size: 0.9em; }
    .pill { display: inline-block; padding: 2px 8px; border-radius: 999px; background: #eef2ff; margin-right: 6px; margin-bottom: 6px; }
    .pathline { display: flex; gap: 8px; align-items: center; }
    .pathline input { flex: 1; }
    .pathline button { width: auto; white-space: nowrap; }
    .browse-list { max-height: 260px; overflow: auto; border: 1px solid #e5e7eb; border-radius: 8px; padding: 8px; background: #fafafa; }
    .browse-item { display: flex; gap: 8px; align-items: center; margin-bottom: 6px; }
    .browse-item button { width: auto; margin: 0; }
  </style>
</head>
<body>
  <h1>MONSDA Web Configurator</h1>
  <p class=\"muted\">Interactive builder for MONSDA configs with samplesheet support.</p>

  <div class=\"card\">
    <h3>1) Samplesheet → SETTINGS</h3>
    <label>Samplesheet path (CSV/TSV)</label>
    <div class="pathline">
      <input id="samplesheetPath" placeholder="/abs/path/to/samplesheet.csv" />
      <button onclick="openPathBrowser('samplesheetPath','samplesheets')">Browse…</button>
    </div>
    <button onclick="parseSamplesheet()">Parse Samplesheet</button>
    <div id="samplesheetStatus"></div>
    <div id="conditionsPreview"></div>
  </div>

  <div class="card">
    <h3>2) Workflow + Tool selection</h3>
    <button onclick="loadTemplate()">Load Workflows from Template</button>
    <div id="templateStatus" class="muted"></div>
    <div id=\"workflowChooser\"></div>
  </div>

  <div class=\"card\">
    <h3>3) Build config</h3>
    <div class=\"row\">
      <div>
        <label>Config name</label>
        <input id=\"configName\" value=\"monsda\" />
      </div>
      <div>
        <label>Output directory</label>
        <div class="pathline">
          <input id=\"outputDir\" placeholder=\"/abs/output/dir\" />
          <button onclick="openPathBrowser('outputDir','dirs')">Browse…</button>
        </div>
      </div>
    </div>
    <label>MAXTHREADS</label>
    <input id=\"maxthreads\" value=\"16\" />

    <label>SETTINGS JSON (auto-filled from samplesheet, editable)</label>
    <textarea id=\"settingsJson\"></textarea>

    <button onclick=\"previewConfig()\">Preview Config JSON</button>
    <button onclick=\"saveConfig()\">Save Config</button>
    <div id=\"buildStatus\"></div>

    <label>Config preview</label>
    <textarea id=\"configPreview\"></textarea>
  </div>

  <div class=\"card\">
    <h3>4) Create project skeleton</h3>
    <label>Project directory</label>
    <div class="pathline">
      <input id="projectDir" placeholder="/abs/path/to/project" />
      <button onclick="openPathBrowser('projectDir','dirs')">Browse…</button>
    </div>
    <button onclick="createProject()">Create Project Structure</button>
    <div id="projectStatus"></div>
  </div>

  <div class="card">
    <h3>Path browser</h3>
    <div id="browserInfo" class="muted"></div>
    <div id="browseRoots" class="muted"></div>
    <div class="pathline">
      <input id="browsePathInput" placeholder="/path/to/browse" />
      <button onclick="browseTo(document.getElementById('browsePathInput').value)">Go</button>
      <button onclick="browseUp()">Up</button>
      <button onclick="chooseCurrentPath()">Use current</button>
    </div>
    <div id="browseStatus" class="muted"></div>
    <div id="browseList" class="browse-list"></div>
  </div>

<script>
let templateFields = null;
let lastConfig = null;
let browserState = { target: '', mode: 'dirs', current: '' };

function setStatus(id, text, ok=true) {
  const el = document.getElementById(id);
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

async function loadTemplate() {
  const holder = document.getElementById('workflowChooser');
  const status = document.getElementById('templateStatus');
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
      const tools = (data[wf] && data[wf].TOOLS) ? Object.keys(data[wf].TOOLS) : [];
      const toolChecks = tools.map(t => `<label><input type=\"checkbox\" data-wf=\"${wf}\" data-tool=\"${t}\" checked /> ${escHtml(t)}</label>`).join('<br>');
      holder.innerHTML += `
        <div style=\"border:1px solid #e5e7eb; border-radius:8px; padding:10px; margin-bottom:8px;\">
          <label><input type=\"checkbox\" data-workflow=\"${wf}\" /> <b>${escHtml(wf)}</b></label>
          <div style=\"margin-left:20px; margin-top:6px;\">${toolChecks || '<span class=\"muted\">No tool list</span>'}</div>
        </div>
      `;
    });

    status.textContent = `Loaded ${workflows.length} workflows from template.`;
    status.className = 'ok';
  } catch (e) {
    holder.innerHTML = '';
    status.textContent = `Failed to load template fields: ${e.message || e}`;
    status.className = 'err';
  }
}

async function loadBrowseRoots() {
  const r = await fetch('/fs/roots');
  const data = await r.json();
  if (!r.ok) throw new Error(data.detail || 'Failed to load roots');
  const rootHtml = (data.roots || []).map(p => `<button onclick=\"browseTo(${JSON.stringify(p)})\">${escHtml(p)}</button>`).join(' ');
  document.getElementById('browseRoots').innerHTML = rootHtml;
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
    document.getElementById('browsePathInput').value = data.current;
    document.getElementById('browserInfo').textContent = `Selecting for: ${browserState.target} (${browserState.mode})`;
    setStatus('browseStatus', `Showing: ${data.current}`, true);

    const dirs = (data.dirs || []).map(d => `
      <div class=\"browse-item\">
        <button onclick=\"browseTo(${JSON.stringify(d.path)})\">📁 ${escHtml(d.name)}</button>
        <button onclick=\"choosePath(${JSON.stringify(d.path)})\">Use</button>
      </div>`).join('');

    const files = (data.files || []).map(f => `
      <div class=\"browse-item\">
        <button onclick=\"choosePath(${JSON.stringify(f.path)})\">📄 ${escHtml(f.name)}</button>
      </div>`).join('');

    document.getElementById('browseList').innerHTML = (dirs + files) || '<div class=\"muted\">No entries.</div>';
  } catch (e) {
    setStatus('browseStatus', e.message || String(e), false);
  }
}

function openPathBrowser(target, mode='dirs') {
  browserState.target = target;
  browserState.mode = mode;
  const start = document.getElementById(target)?.value?.trim() || '';
  loadBrowseRoots().then(() => browseTo(start));
}

function choosePath(path) {
  if (!browserState.target) return;
  document.getElementById(browserState.target).value = path;
  setStatus('browseStatus', `Selected: ${path}`, true);
}

function chooseCurrentPath() {
  if (!browserState.current) return;
  choosePath(browserState.current);
}

function browseUp() {
  if (!browserState.current) return;
  const parent = browserState.current === '/' ? '/' : browserState.current.replace(/\/+$/, '').replace(/\/[^\/]*$/, '') || '/';
  browseTo(parent);
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
      document.getElementById('conditionsPreview').innerHTML = data.conditions.map(c => `<span class=\"pill\">${c}</span>`).join('');
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
    tools[wf] = Array.from(document.querySelectorAll(`input[data-wf="${wf}"][data-tool]:checked`)).map(i => i.getAttribute('data-tool'));
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
    const workflows = (lastConfig && lastConfig.WORKFLOWS) ? lastConfig.WORKFLOWS.split(',').map(x => x.trim()).filter(Boolean) : [];
    const r = await fetch('/project/create', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        project_dir: document.getElementById('projectDir').value.trim(),
        workflows: workflows
      })
    });
    const data = await r.json();
    if (r.ok) {
      setStatus('projectStatus', `Created: ${data.path}`);
    } else {
      setStatus('projectStatus', data.detail || 'Failed to create project.', false);
    }
  } catch (e) {
    setStatus('projectStatus', 'Failed to create project.', false);
  }
}

window.onload = loadTemplate;
</script>
</body>
</html>
"""
