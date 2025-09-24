import json
import os
from typing import Any, Dict

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from pydantic import BaseModel

TEMPLATE_PATH = os.path.join(
    os.path.dirname(__file__), "../configs/template_base_commented.json"
)

app = FastAPI(title="MONSDA Configurator Web Service")

# Allow CORS for local development
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


class ConfigRequest(BaseModel):
    config_name: str
    config: Dict[str, Any]
    output_dir: str


def load_template() -> Dict[str, Any]:
    with open(TEMPLATE_PATH, "r") as f:
        return json.load(f)


def strip_comments(d):
    if isinstance(d, dict):
        return {k: strip_comments(v) for k, v in d.items() if k != "comment"}
    elif isinstance(d, list):
        return [strip_comments(x) for x in d]
    else:
        return d

    config_name: str
    config: Dict[str, Any]
    output_dir: str


@app.get("/template", response_model=Dict[str, Any])
def get_template():
    """Get the config template (with comments)."""
    return load_template()


@app.get("/template/fields", response_model=Dict[str, Any])
def get_template_fields():
    """Get the config template (without comments)."""
    return strip_comments(load_template())


@app.post("/generate_config")
def generate_config(req: ConfigRequest):
    """Generate a config file from user input. User specifies output_dir."""
    config_name = req.config_name.strip()
    output_dir = os.path.abspath(req.output_dir.strip())
    if not config_name or any(c in config_name for c in "/\\"):
        raise HTTPException(status_code=400, detail="Invalid config name.")
    if not output_dir or not os.path.isdir(output_dir):
        raise HTTPException(status_code=400, detail="Invalid output directory.")
    config_path = os.path.join(output_dir, f"config_{config_name}.json")
    with open(config_path, "w") as f:
        json.dump(req.config, f, indent=4)
    return {"message": "Config generated.", "path": config_path}


# Download config by full path (for security, restrict to files under allowed parent dir)
@app.get("/download_config/")
def download_config(config_path: str):
    config_path = os.path.abspath(config_path)
    if not os.path.exists(config_path):
        raise HTTPException(status_code=404, detail="Config not found.")
    if not config_path.endswith(".json"):
        raise HTTPException(status_code=400, detail="Only .json files allowed.")
    return FileResponse(config_path, filename=os.path.basename(config_path))


class DirRequest(BaseModel):
    config: Dict[str, Any]
    project_dir: str


def safe_makedirs(path):
    os.makedirs(path, exist_ok=True)


def create_project_structure(config: Dict[str, Any], project_dir: str):
    # Example: create folders for workflows, logs, counts, etc.
    safe_makedirs(project_dir)
    for wf in config.get("WORKFLOWS", "").split(","):
        wf = wf.strip()
        if wf:
            safe_makedirs(os.path.join(project_dir, wf))
    safe_makedirs(os.path.join(project_dir, "LOGS"))
    safe_makedirs(os.path.join(project_dir, "COUNTS"))
    # Add more as needed


@app.post("/generate_project_dir")
def generate_project_dir(req: DirRequest):
    project_dir = os.path.abspath(req.project_dir)
    if not project_dir.startswith(os.getcwd()):
        raise HTTPException(
            status_code=400, detail="Project dir must be inside working directory."
        )
    if os.path.exists(project_dir) and os.listdir(project_dir):
        raise HTTPException(
            status_code=400, detail="Project dir already exists and is not empty."
        )
    create_project_structure(req.config, project_dir)
    return {"message": "Project directory structure created.", "path": project_dir}


from fastapi.responses import HTMLResponse


@app.get("/", response_class=HTMLResponse)
def root():
    return """
<!DOCTYPE html>
<html lang=\"en\">
<head>
    <meta charset=\"UTF-8\">
    <title>MONSDA Configurator</title>
    <style>
        body { font-family: sans-serif; margin: 2em; background: #f8f8fa; }
        h1 { color: #2c3e50; }
        textarea, input, select { width: 100%; margin: 0.5em 0; }
        .section { background: #fff; border-radius: 8px; box-shadow: 0 2px 8px #0001; padding: 1.5em; margin-bottom: 2em; }
        button { background: #2c3e50; color: #fff; border: none; padding: 0.7em 1.5em; border-radius: 4px; cursor: pointer; }
        button:disabled { background: #aaa; }
        .success { color: green; }
        .error { color: red; }
        .field-label { font-weight: bold; margin-top: 1em; }
    </style>
</head>
<body>
    <h1>MONSDA Configurator</h1>
    <div class=\"section\">
        <h2>Step 1: Edit Configuration</h2>
        <button onclick=\"loadTemplate()\">Load Template</button>
        <span id=\"template-status\"></span>
        <textarea id=\"config-editor\" rows=\"24\" placeholder=\"Config JSON will appear here...\"></textarea>
    </div>
    <div class=\"section\">
        <h2>Step 2: Generate Config File</h2>
        <label class=\"field-label\">Config Name (no spaces or slashes):</label>
        <input id=\"config-name\" type=\"text\" placeholder=\"e.g. testproject\" />
        <label class=\"field-label\">Config Output Directory (absolute path):</label>
        <input id=\"config-output-dir\" type=\"text\" placeholder=\"e.g. /home/user/configs\" />
        <button onclick=\"generateConfig()\">Generate Config</button>
        <span id=\"generate-status\"></span>
        <div id=\"download-link\"></div>
    </div>
    <div class=\"section\">
        <h2>Step 3: Create Project Directory</h2>
        <label class=\"field-label\">Project Directory (absolute path):</label>
        <input id=\"project-dir\" type=\"text\" placeholder=\"e.g. /home/user/myproject\" />
        <button onclick=\"createProjectDir()\">Create Directory Structure</button>
        <span id=\"dir-status\"></span>
    </div>
    <script>
    function loadTemplate() {
        fetch('/template/fields').then(r => r.json()).then(data => {
            document.getElementById('config-editor').value = JSON.stringify(data, null, 4);
            document.getElementById('template-status').textContent = 'Template loaded.';
        }).catch(e => {
            document.getElementById('template-status').textContent = 'Failed to load template.';
        });
    }

    function generateConfig() {
        const configText = document.getElementById('config-editor').value;
        const configName = document.getElementById('config-name').value.trim();
        const outputDir = document.getElementById('config-output-dir').value.trim();
        let config;
        try {
            config = JSON.parse(configText);
        } catch (e) {
            document.getElementById('generate-status').textContent = 'Invalid JSON.';
            document.getElementById('generate-status').className = 'error';
            return;
        }
        fetch('/generate_config', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ config_name: configName, config: config, output_dir: outputDir })
        }).then(r => r.json()).then(data => {
            if (data.path) {
                document.getElementById('generate-status').textContent = 'Config generated!';
                document.getElementById('generate-status').className = 'success';
                document.getElementById('download-link').innerHTML = `<a href='/download_config/?config_path=${encodeURIComponent(data.path)}' target='_blank'>Download ${data.path.split('/').pop()}</a>`;
            } else {
                document.getElementById('generate-status').textContent = data.detail || 'Error generating config.';
                document.getElementById('generate-status').className = 'error';
            }
        }).catch(e => {
            document.getElementById('generate-status').textContent = 'Error generating config.';
            document.getElementById('generate-status').className = 'error';
        });
    }

    function createProjectDir() {
        const configText = document.getElementById('config-editor').value;
        const projectDir = document.getElementById('project-dir').value.trim();
        let config;
        try {
            config = JSON.parse(configText);
        } catch (e) {
            document.getElementById('dir-status').textContent = 'Invalid JSON.';
            document.getElementById('dir-status').className = 'error';
            return;
        }
        fetch('/generate_project_dir', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ config: config, project_dir: projectDir })
        }).then(r => r.json()).then(data => {
            if (data.path) {
                document.getElementById('dir-status').textContent = 'Project directory created!';
                document.getElementById('dir-status').className = 'success';
            } else {
                document.getElementById('dir-status').textContent = data.detail || 'Error creating directory.';
                document.getElementById('dir-status').className = 'error';
            }
        }).catch(e => {
            document.getElementById('dir-status').textContent = 'Error creating directory.';
            document.getElementById('dir-status').className = 'error';
        });
    }
    // Auto-load template on page load
    window.onload = loadTemplate;
    </script>
</body>
</html>
"""
