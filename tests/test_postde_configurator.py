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
TEMPLATE = REPO / "configs" / "template_base_commented.json"

# --- web configurator ------------------------------------------------------

try:
    from fastapi.testclient import TestClient

    from MONSDA import web_configurator as wc
    from MONSDA.extras.PostDE import POSTDE_DEFAULTS, load_postde_config

    HAS_HTTPX = True
except ImportError:  # pragma: no cover - fallback for envs without httpx
    HAS_HTTPX = False
    from MONSDA import web_configurator as wc
    from MONSDA.extras.PostDE import POSTDE_DEFAULTS, load_postde_config


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
    return {"enabled": True}


def _client():
    return TestClient(wc.app)


@pytest.mark.skipif(not HAS_HTTPX, reason="httpx not available")
class TestWebPreview:
    def test_preview_attaches_disabled_postde(self):
        r = _client().post("/config/preview", json=_preview_body())
        assert r.status_code == 200
        assert r.json()["config"]["DE"]["POSTDE"] == {"enabled": False}
        assert "POSTDE" not in r.json()["config"]

    def test_preview_enabled_postde_has_defaults(self):
        r = _client().post(
            "/config/preview", json=_preview_body(postde=_enabled_postde())
        )
        assert r.status_code == 200
        section = r.json()["config"]["DE"]["POSTDE"]
        assert section == POSTDE_DEFAULTS
        assert section["background"]["expressed"]["metric"] == "cpm"

    def test_preview_without_de_has_no_postde(self):
        r = _client().post(
            "/config/preview",
            json=_preview_body(workflows=("QC",), postde=_enabled_postde()),
        )
        assert r.status_code == 200
        assert "POSTDE" not in r.json()["config"]
        assert "POSTDE" not in r.json()["config"].get("QC", {})


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

    def _save(self, tmp_path, cfg):
        return _client().post(
            "/config/save",
            json={"config_name": "t", "output_dir": str(tmp_path), "config": cfg},
        )

    def test_save_validates_postde_before_write(self, tmp_path):
        cfg = self._valid_cfg(tmp_path)
        cfg["DE"]["POSTDE"]["background"] = {"file": "missing.tsv"}
        r = self._save(tmp_path, cfg)
        assert r.status_code == 400
        assert not list(tmp_path.glob("config_*.json"))

    def test_save_rejects_legacy_config_key(self, tmp_path):
        cfg = self._valid_cfg(tmp_path)
        cfg["DE"]["POSTDE"]["config"] = "old.json"
        r = self._save(tmp_path, cfg)
        assert r.status_code == 400
        assert "unknown top-level" in r.json()["detail"]

    def test_save_unsupported_de_tool_raises(self, tmp_path):
        cfg = self._valid_cfg(tmp_path)
        cfg["DE"]["TOOLS"] = {"diego": "Analysis/DE/DIEGO.R"}
        r = self._save(tmp_path, cfg)
        assert r.status_code == 400
        assert "does not support" in r.json()["detail"]

    def test_save_writes_postde(self, tmp_path):
        cfg = self._valid_cfg(tmp_path)
        r = self._save(tmp_path, cfg)
        assert r.status_code == 200
        written = json.loads((tmp_path / "config_t.json").read_text())
        assert written["DE"]["POSTDE"] == POSTDE_DEFAULTS

    def test_create_validates_postde_before_creation(self, tmp_path):
        cfg = self._valid_cfg(tmp_path)
        cfg["DE"]["POSTDE"]["background"] = {"file": "missing.tsv"}
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
        cfg = self._valid_cfg(tmp_path)
        path = tmp_path / "config_t.json"
        path.write_text(json.dumps(cfg))
        reloaded = json.loads(path.read_text())
        section = load_postde_config(reloaded)
        assert section["enabled"] is True
        assert section["analysis"]["background"] == POSTDE_DEFAULTS["background"]


@pytest.mark.skipif(not HAS_HTTPX, reason="httpx not available")
class TestWebHtml:
    def test_served_html_has_postde_controls(self):
        html = _client().get("/").text
        assert "postdeEnabled" in html
        assert "postdeConfig" not in html

    def test_served_html_does_not_offer_postde_as_workflow(self):
        html = _client().get("/").text
        assert "data-workflow=\"POSTDE\"" not in html


class TestTemplate:
    def test_template_has_disabled_de_postde(self):
        template = wc.strip_comments(wc.load_template())
        assert "POSTDE" not in template
        assert template["DE"]["POSTDE"]["enabled"] is False

    def test_template_postde_validates_when_enabled(self):
        template = wc.strip_comments(wc.load_template())
        template["WORKFLOWS"] = "DE"
        template["DE"]["POSTDE"]["enabled"] = True
        assert load_postde_config(template)["enabled"] is True


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
    def test_no_de_leaves_config_untouched(self, monkeypatch):
        project, _ = _stub_guide(monkeypatch, [])
        final_dict = {"WORKFLOWS": "QC"}
        configurator.configure_postde(final_dict)
        assert final_dict == {"WORKFLOWS": "QC"}

    def test_choice_no_disables(self, monkeypatch):
        _stub_guide(monkeypatch, ["no"])
        final_dict = {"WORKFLOWS": "DE", "DE": {}}
        configurator.configure_postde(final_dict)
        assert final_dict["DE"]["POSTDE"] == {"enabled": False}

    def test_choice_yes_writes_defaults(self, monkeypatch):
        _stub_guide(monkeypatch, ["yes"])
        final_dict = {
            "WORKFLOWS": "DE",
            "DE": {"TOOLS": {"deseq2": "Analysis/DE/DESEQ2.R"}},
        }
        configurator.configure_postde(final_dict)
        assert final_dict["DE"]["POSTDE"] == POSTDE_DEFAULTS
        assert load_postde_config(final_dict)["enabled"] is True

    def test_unsupported_de_tool_reported(self, monkeypatch, capsys):
        _stub_guide(monkeypatch, ["yes"])
        final_dict = {
            "WORKFLOWS": "DE",
            "DE": {"TOOLS": {"diego": "Analysis/DE/DIEGO.R"}},
        }
        configurator.configure_postde(final_dict)
        assert "does not support DE tools" in capsys.readouterr().out

    def test_old_pickled_project_without_postde_attr(self, monkeypatch):
        project, _ = _stub_guide(monkeypatch, ["no"])
        del project.postde
        final_dict = {"WORKFLOWS": "DE", "DE": {}}
        configurator.configure_postde(final_dict)
        assert final_dict["DE"]["POSTDE"] == {"enabled": False}

    def _old_config(self, tmp_path, postde):
        de = {"TOOLS": {"deseq2": "Analysis/DE/DESEQ2.R"}}
        if postde is not None:
            de["POSTDE"] = postde
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
                    "DE": de,
                }
            )
        )
        return cfg_path

    def test_modify_config_without_postde(self, tmp_path, monkeypatch):
        cfg_path = self._old_config(tmp_path, None)
        project, _ = _stub_guide(monkeypatch, ["5", "no", ""])
        project.mode = "modify"
        captured = {}
        monkeypatch.setattr(
            configurator, "create_project", lambda d: captured.update(final_dict=d)
        )
        configurator.modify(str(cfg_path))
        assert captured["final_dict"]["DE"]["POSTDE"] == {"enabled": False}

    def test_modify_keeps_existing_postde(self, tmp_path, monkeypatch):
        section = {"enabled": True, "padj": 0.01, "gprofiler": {"enabled": False}}
        cfg_path = self._old_config(tmp_path, section)
        project, _ = _stub_guide(monkeypatch, ["5", "yes", ""])
        project.mode = "modify"
        captured = {}
        monkeypatch.setattr(
            configurator, "create_project", lambda d: captured.update(final_dict=d)
        )
        configurator.modify(str(cfg_path))
        out = captured["final_dict"]["DE"]["POSTDE"]
        assert out["enabled"] is True
        assert out["padj"] == 0.01
