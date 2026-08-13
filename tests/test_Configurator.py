import importlib
import sys
from pathlib import Path

import pytest

pytest.importorskip("snakemake")

import snakemake.common.configfile as configfile

REPO = Path(__file__).resolve().parents[1]
ORIGINAL_LOAD_CONFIGFILE = configfile.load_configfile


def _load_repo_template(_):
    return ORIGINAL_LOAD_CONFIGFILE(
        str(REPO / "configs" / "template_base_commented.json")
    )


configfile.load_configfile = _load_repo_template
original_argv = sys.argv
sys.argv = [sys.argv[0]]
configurator = importlib.import_module("MONSDA.Configurator")
sys.argv = original_argv
configfile.load_configfile = ORIGINAL_LOAD_CONFIGFILE


@pytest.mark.parametrize(
    ("tool_answer", "expect_features"),
    [("1", True), ("2", False), ("1,2", True)],
)
def test_counting_features_only_prompted_for_countreads(
    monkeypatch, tool_answer, expect_features
):
    project = configurator.PROJECT()
    project.baseDict["COUNTING"] = {
        "TOOLS": {
            "countreads": "featureCounts",
            "salmon": "salmon",
            "oarfish": "oarfish",
            "alevinfry": "alevin-fry",
            "simpleaf": "simpleaf",
        },
        "FEATURES": {"exon": "gene_id", "gene": "gene_id"},
    }
    project.workflowsDict["COUNTING"]
    guide = configurator.GUIDE()
    prompts = []

    def display(options=None, question=None, **_):
        prompts.append((question, dict(options or {})))
        if question.startswith("Select from these available Tools"):
            configurator.guide.answer = tool_answer
        elif question.startswith("Select from these available FEATURES"):
            configurator.guide.answer = "1,2"

    monkeypatch.setattr(configurator, "project", project, raising=False)
    monkeypatch.setattr(configurator, "guide", guide, raising=False)
    monkeypatch.setattr(configurator, "pickle_unfinished", lambda _: None)
    monkeypatch.setattr(configurator, "show_settings", lambda: None)
    monkeypatch.setattr(configurator.guide, "display", display)
    monkeypatch.setattr(configurator.guide, "clear", lambda *_: None)

    configurator.set_workflows("COUNTING")

    feature_prompts = [entry for entry in prompts if "FEATURES" in entry[0]]
    assert bool(feature_prompts) is expect_features
    if expect_features:
        assert feature_prompts[0][1] == {1: "exon", 2: "gene"}
