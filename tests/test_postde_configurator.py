"""Focused tests for the PostDE slice of the web configurator and the CLI.

Covers the web configurator (preview/save/create validation, disabled
backwards compatibility, served HTML, serialized config reload) and the CLI
``configure_postde`` function (choices and DE-tool exclusion). Uses tmp_path
only, no local hardcoded paths.
"""

import importlib
import json
import os
import sys
from pathlib import Path

import pytest

pytest.importorskip("snakemake")

REPO = Path(__file__).resolve().parents[1]
EXAMPLE = REPO / "configs" / "postde_analysis.json"
TEMPLATE = REPO / "configs" / "template_base_commented.json"

# --- web configurator ------------------------------------------------------

try:
    from fastapi.testclient import TestClient

    from MONSDA import web_configurator as wc
    from MONSDA.PostDE import load_postde_config

    HAS_HTTPX = True
except ImportError:  # pragma: no cover - fallback for envs without httpx
    HAS_HTTPX = False
    from MONSDA import web_configurator as wc
    from MONSDA.PostDE import load_postde_config


def _settings():
    return {"Exp": {"WT": {"SAMPLES": ["s1", "s2"], "GROUPS": ["g", "g"]}}}


def _preview_body(workflows=("DE",), postde=None, tools=None):
    body = {
        "config_name": "t",
        "output_dir": str(REPO),
        "workflows": list(workflows),
        "tools": tools or {},
        "maxthreads": "16",
        "settings": _settings(),
        "samplesheet_path": None,
    }
    if postde is not None:
        body["postde"] = postde
    return body


def _enabled_postde():
    return {"enabled": True, "config": str(EXAMPLE)}


def _client():
    return TestClient(wc.app)


@pytest.mark.skipif(not HAS_HTTPX, reason="httpx not available")
class TestWebPreview:
    def test_preview_attaches_disabled_postde(self):
        r = _client().post("/config/preview", json=_preview_body())
        assert r.status_code == 200
        assert r.json()["config"]["POSTDE"] == {"enabled": False, "config": ""}

    def test_preview_enabled_postde_validates(self):
        r = _client().post(
            "/config/preview", json=_preview_body(postde=_enabled_postde())
        )
        assert r.status_code == 200
        section = r.json()["config"]["POSTDE"]
        assert section["enabled"] is True
        assert section["config"] == str(EXAMPLE.resolve())

    def test_preview_disabled_postde_backwards_compat(self):
        r = _client().post(
            "/config/preview",
            json=_preview_body(postde={"enabled": False, "config": ""}),
        )
        assert r.status_code == 200
        assert r.json()["config"]["POSTDE"] == {"enabled": False, "config": ""}

    def test_preview_missing_de_raises(self):
        r = _client().post(
            "/config/preview",
            json=_preview_body(workflows=("QC",), postde=_enabled_postde()),
        )
        assert r.status_code == 400
        assert "DE in WORKFLOWS" in r.json()["detail"]

    def test_preview_invalid_postde_section_raises(self):
        r = _client().post(
            "/config/preview",
            json=_preview_body(postde={"enabled": "yes", "config": ""}),
        )
        assert r.status_code == 400


@pytest.mark.skipif(not HAS_HTTPX, reason="httpx not available")
class TestWebSaveCreate:
    def _valid_cfg(self, tmp_path):
        return wc.build_config(
            wc.BuildConfigRequest(
                config_name="t",
                output_dir=str(tmp_path),
                workflows=["DE"],
                settings=_settings(),
                postde=_enabled_postde(),
            )
        )

    def test_save_validates_postde_before_write(self, tmp_path):
        cfg = self._valid_cfg(tmp_path)
        cfg["POSTDE"] = {"enabled": True, "config": "missing.json"}
        r = _client().post(
            "/config/save",
            json={
                "config_name": "t",
                "output_dir": str(tmp_path),
                "config": cfg,
            },
        )
        assert r.status_code == 400
        assert not list(tmp_path.glob("config_*.json"))

    def test_save_unsupported_de_tool_raises(self, tmp_path):
        cfg = self._valid_cfg(tmp_path)
        cfg["DE"]["TOOLS"] = {"diego": "Analysis/DE/DIEGO.R"}
        r = _client().post(
            "/config/save",
            json={
                "config_name": "t",
                "output_dir": str(tmp_path),
                "config": cfg,
            },
        )
        assert r.status_code == 400
        assert "does not support" in r.json()["detail"]

    def test_save_writes_normalized_postde(self, tmp_path):
        cfg = self._valid_cfg(tmp_path)
        r = _client().post(
            "/config/save",
            json={
                "config_name": "t",
                "output_dir": str(tmp_path),
                "config": cfg,
            },
        )
        assert r.status_code == 200
        written = json.loads((tmp_path / "config_t.json").read_text())
        assert written["POSTDE"]["enabled"] is True
        assert written["POSTDE"]["config"] == str(EXAMPLE.resolve())

    def test_create_validates_postde_before_creation(self, tmp_path):
        cfg = self._valid_cfg(tmp_path)
        cfg["POSTDE"] = {"enabled": True, "config": "missing.json"}
        r = _client().post(
            "/project/create",
            json={
                "project_dir": str(tmp_path),
                "project_name": "proj",
                "condition_files": [],
                "config": cfg,
                "settings": _settings(),
            },
        )
        assert r.status_code == 400
        assert not (tmp_path / "proj").exists()

    def test_create_unsupported_de_tool_raises(self, tmp_path):
        cfg = self._valid_cfg(tmp_path)
        cfg["DE"]["TOOLS"] = {"diego": "Analysis/DE/DIEGO.R"}
        r = _client().post(
            "/project/create",
            json={
                "project_dir": str(tmp_path),
                "project_name": "proj",
                "condition_files": [],
                "config": cfg,
                "settings": _settings(),
            },
        )
        assert r.status_code == 400
        assert not (tmp_path / "proj").exists()

    def test_serialized_config_reloads(self, tmp_path):
        cfg = wc.build_config(
            wc.BuildConfigRequest(
                config_name="t",
                output_dir=str(tmp_path),
                workflows=["DE"],
                settings=_settings(),
                postde=_enabled_postde(),
            )
        )
        path = tmp_path / "config_t.json"
        path.write_text(json.dumps(cfg))
        reloaded = json.loads(path.read_text())
        section = load_postde_config(reloaded["POSTDE"], config=reloaded)
        assert section["enabled"] is True
        assert section["config"] == str(EXAMPLE.resolve())


@pytest.mark.skipif(not HAS_HTTPX, reason="httpx not available")
class TestWebHtml:
    def test_served_html_has_postde_controls(self):
        html = _client().get("/").text
        assert "postdeEnabled" in html
        assert "postdeConfig" in html
        assert "openPathBrowser('postdeConfig','all')" in html

    def test_served_html_does_not_offer_postde_as_workflow(self):
        html = _client().get("/").text
        # the frontend workflow filter must exclude POSTDE from the chooser
        assert "'POSTDE'" in html
        assert "data-workflow=\"POSTDE\"" not in html


class TestTemplateAndExample:
    def test_template_has_disabled_postde(self):
        template = wc.strip_comments(wc.load_template())
        assert template["POSTDE"] == {"enabled": False, "config": ""}

    def test_example_validates_from_arbitrary_cwd(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        section = load_postde_config({"enabled": True, "config": str(EXAMPLE)})
        assert section["enabled"] is True
        assert section["config"] == str(EXAMPLE.resolve())

    def test_example_requires_no_external_resources(self):
        analysis = json.loads(EXAMPLE.read_text())
        for method, keys in (
            ("gprofiler", ("background",)),
            ("clusterprofiler", ("term2gene", "term2name", "background")),
            ("decoupler", ("network",)),
            ("dream", ("metadata",)),
        ):
            for key in keys:
                assert analysis[method].get(key) in (None, ""), (
                    method + "." + key + " must not require a resource file"
                )


# --- CLI configurator ------------------------------------------------------

import snakemake.common.configfile as configfile  # noqa: E402

_ORIGINAL_LOAD_CONFIGFILE = configfile.load_configfile


def _load_repo_template(_):
    return _ORIGINAL_LOAD_CONFIGFILE(str(TEMPLATE))


configfile.load_configfile = _load_repo_template
_original_argv = sys.argv
sys.argv = [sys.argv[0]]
configurator = importlib.import_module("MONSDA.Configurator")
sys.argv = _original_argv
configfile.load_configfile = _ORIGINAL_LOAD_CONFIGFILE
# the import-time monkeypatch above bound Configurator's load_configfile to the
# template loader; restore it so modify() reads real config files
configurator.load_configfile = _ORIGINAL_LOAD_CONFIGFILE


def _stub_guide(monkeypatch, answers):
    project = configurator.PROJECT()
    guide = configurator.GUIDE()
    monkeypatch.setattr(configurator, "project", project, raising=False)
    monkeypatch.setattr(configurator, "guide", guide, raising=False)
    monkeypatch.setattr(configurator, "pickle_unfinished", lambda _: None)
    it = iter(answers)

    def fake_display(options=None, question=None, proof=None, spec=None, whitespace=False):
        configurator.guide.answer = next(it)

    monkeypatch.setattr(configurator.guide, "display", fake_display)
    monkeypatch.setattr(configurator.guide, "clear", lambda *_: None)
    return project, guide


class TestCliConfigurePostde:
    def test_no_de_outputs_disabled(self, monkeypatch):
        project, _ = _stub_guide(monkeypatch, [])
        project.workflowsDict["QC"]
        final_dict = {"WORKFLOWS": "QC"}
        configurator.configure_postde(final_dict)
        assert final_dict["POSTDE"] == {"enabled": False, "config": ""}

    def test_choice_no_disables(self, monkeypatch):
        project, _ = _stub_guide(monkeypatch, ["no"])
        project.workflowsDict["DE"]
        final_dict = {"WORKFLOWS": "DE"}
        configurator.configure_postde(final_dict)
        assert final_dict["POSTDE"] == {"enabled": False, "config": ""}

    def test_choice_yes_enables_with_abs_path(self, monkeypatch):
        project, _ = _stub_guide(monkeypatch, ["yes", str(EXAMPLE)])
        project.workflowsDict["DE"]
        final_dict = {"WORKFLOWS": "DE"}
        configurator.configure_postde(final_dict)
        assert final_dict["POSTDE"]["enabled"] is True
        assert final_dict["POSTDE"]["config"] == str(EXAMPLE.resolve())

    def test_invalid_path_retries(self, monkeypatch):
        project, _ = _stub_guide(monkeypatch, ["yes", "missing.json", str(EXAMPLE)])
        project.workflowsDict["DE"]
        final_dict = {"WORKFLOWS": "DE"}
        configurator.configure_postde(final_dict)
        assert final_dict["POSTDE"]["enabled"] is True

    def test_unsupported_de_tool_rejected(self, monkeypatch, capsys):
        project, _ = _stub_guide(monkeypatch, ["yes", str(EXAMPLE)])
        project.workflowsDict["DE"]
        final_dict = {
            "WORKFLOWS": "DE",
            "DE": {"TOOLS": {"diego": "Analysis/DE/DIEGO.R"}},
        }
        # the error is reported and the path prompt is retried
        with pytest.raises(StopIteration):
            configurator.configure_postde(final_dict)
        assert "does not support DE tools" in capsys.readouterr().out

    def test_old_pickled_project_without_postde_attr(self, monkeypatch):
        project, _ = _stub_guide(monkeypatch, ["no"])
        project.workflowsDict["DE"]
        del project.postde  # simulate a pickled project from before PostDE
        final_dict = {"WORKFLOWS": "DE"}
        configurator.configure_postde(final_dict)
        assert final_dict["POSTDE"] == {"enabled": False, "config": ""}

    def test_modify_old_config_without_postde_key(self, tmp_path, monkeypatch):
        cfg_path = tmp_path / "config_old.json"
        cfg_path.write_text(
            json.dumps(
                {
                    "WORKFLOWS": "DE",
                    "MAXTHREADS": "8",
                    "VERSION": "x",
                    "BINS": "MONSDA/scripts",
                    "SETTINGS": {
                        "Exp": {"WT": {"SAMPLES": ["s1"], "GROUPS": ["g"]}}
                    },
                    "DE": {"TOOLS": {"deseq2": "Analysis/DE/DESEQ2.R"}},
                }
            )
        )
        project, _ = _stub_guide(monkeypatch, ["5", "no", ""])
        project.mode = "modify"
        captured = {}

        def fake_create_project(final_dict):
            captured["final_dict"] = final_dict

        monkeypatch.setattr(configurator, "create_project", fake_create_project)
        configurator.modify(str(cfg_path))
        assert captured["final_dict"]["POSTDE"] == {"enabled": False, "config": ""}
        assert "DE" in captured["final_dict"]["WORKFLOWS"]
