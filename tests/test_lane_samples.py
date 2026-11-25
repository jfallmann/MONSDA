import glob
import gzip
import json
import os
import re
import shutil
import subprocess
import types
from pathlib import Path

import pytest

pytest.importorskip("snakemake")

import MONSDA  # noqa: E402
import MONSDA.Params as mp  # noqa: E402
import MONSDA.RunMONSDA as mr  # noqa: E402
import MONSDA.Workflows as mw  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "tests" / "data"

SAMPLES = ["S1", "S2", "S3", "S8", "S9", "S10"]
LANES = ["L001", "L002"]
READS = ["R1", "R2"]
FORM = {
    "S1": {"L001": "A", "L002": "A"},
    "S2": {"L001": "A", "L002": "A"},
    "S3": {"L001": "B", "L002": "A"},
    "S8": {"L001": "B", "L002": "B"},
    "S9": {"L001": "B", "L002": "B"},
    "S10": {"L001": "B", "L002": "B"},
}


def _lane_name(sample, lane, read):
    if FORM[sample][lane] == "A":
        return f"{sample}_{lane}_{read}_001.fastq.gz"
    return f"{sample}_{read}_{lane}_001.fastq.gz"


def _lane_text(sample, lane, read):
    return f"{sample} lane {lane} read {read}\n"


def _write_fastq(directory, name, text):
    directory.mkdir(parents=True, exist_ok=True)
    with gzip.open(directory / name, "wt") as fh:
        fh.write(text)


def _read_fastq(path):
    with gzip.open(path, "rt") as fh:
        return fh.read()


@pytest.fixture
def workdir(tmp_path, monkeypatch):
    fastq = tmp_path / "FASTQ" / "Test"
    fastq.mkdir(parents=True)
    for sample in SAMPLES:
        for lane in LANES:
            for read in READS:
                with gzip.open(fastq / _lane_name(sample, lane, read), "wt") as fh:
                    fh.write(_lane_text(sample, lane, read))
    cfg = json.loads((DATA / "config_Test.json").read_text())
    cfg["MAXTHREADS"] = "4"
    cfg["SETTINGS"]["Test"]["SAMPLES"] = list(SAMPLES)
    cfg["SETTINGS"]["Test"]["SEQUENCING"] = "paired"
    monkeypatch.setattr(mw, "workflowpath", str(REPO / "workflows"))
    monkeypatch.setattr(mw, "envpath", str(REPO / "envs") + os.sep)
    old = os.getcwd()
    os.chdir(tmp_path)
    try:
        yield cfg, fastq
    finally:
        os.chdir(old)


def test_lane_merge_combined_and_reverse_order(workdir):
    cfg, fastq = workdir
    raw = {p.name: p.read_bytes() for p in sorted(fastq.glob("*.fastq.gz"))}
    assert len(raw) == 24

    assert mp.prepare_lane_split_fastqs(cfg) == 12

    for sample in SAMPLES:
        for read in READS:
            target = fastq / f"{sample}_{read}.fastq.gz"
            with gzip.open(target, "rt") as fh:
                assert fh.read() == "".join(
                    _lane_text(sample, lane, read) for lane in LANES
                ), target.name
    assert {p.name: p.read_bytes() for p in fastq.glob("*.fastq.gz") if p.name in raw} == raw

    assert mp.prepare_lane_split_fastqs(cfg) == 0


def _expected_merge(fastq, names, read):
    ordered = sorted(
        (n for n in names if mp._parse_fastq_name(n, "S1").read == read),
        key=lambda n: (mp._split_index(n), n),
    )
    return "".join(_read_fastq(fastq / n) for n in ordered)


@pytest.mark.parametrize(
    "sequencing,names,targets",
    [
        (
            "paired",
            [
                "S1_L002_R1_001.fastq.gz",
                "S1_L001_R1_001.fastq.gz",
                "S1_R2_L002_001.fastq.gz",
                "S1_R2_L001_001.fastq.gz",
            ],
            ["S1_R1.fastq.gz", "S1_R2.fastq.gz"],
        ),
        (
            "singlecell",
            [
                "S1_L001_R1_001.fastq.gz",
                "S1_L002_R1_001.fastq.gz",
                "S1_L001_R2_001.fastq.gz",
                "S1_L002_R2_001.fastq.gz",
            ],
            ["S1_R1.fastq.gz", "S1_R2.fastq.gz"],
        ),
        (
            "single",
            ["S1_L002_R1_001.fastq.gz", "S1_L001_R1_001.fastq.gz"],
            ["S1.fastq.gz"],
        ),
    ],
)
def test_prepare_materializes_canonical_layouts(
    fastq_cwd, sequencing, names, targets
):
    fastq = fastq_cwd / "FASTQ" / "Test"
    for name in names:
        _write_fastq(fastq, name, f"{name}\n")
    raw = {p.name: p.read_bytes() for p in fastq.iterdir()}

    cfg = _settings(["S1"], sequencing=sequencing)
    assert mp.prepare_lane_split_fastqs(cfg) == len(targets)

    for target in targets:
        read = int(re.search(r"_R(\d)", target).group(1)) if "_R" in target else 1
        assert _read_fastq(fastq / target) == _expected_merge(fastq, names, read)
    assert {p.name: p.read_bytes() for p in fastq.iterdir() if p.name in raw} == raw
    assert mp.prepare_lane_split_fastqs(cfg) == 0


def test_prepare_lowercase_mates_and_chunks(fastq_cwd):
    fastq = fastq_cwd / "FASTQ" / "Test"
    _write_fastq(fastq, "S1_r1.fastq.gz", "r1\n")
    _write_fastq(fastq, "S1_R2_001.fastq.gz", "r2 chunk\n")
    cfg = _settings(["S1"])
    assert mp.prepare_lane_split_fastqs(cfg) == 2
    assert _read_fastq(fastq / "S1_R1.fastq.gz") == "r1\n"
    assert _read_fastq(fastq / "S1_R2.fastq.gz") == "r2 chunk\n"
    assert mp.get_samples(cfg) == ["Test/S1_R1", "Test/S1_R2"]
    assert mp.prepare_lane_split_fastqs(cfg) == 0


def test_prepare_existing_canonical_wins(fastq_cwd):
    fastq = fastq_cwd / "FASTQ" / "Test"
    _write_fastq(fastq, "S1_R1.fastq.gz", "canonical R1\n")
    _write_fastq(fastq, "S1_R2.fastq.gz", "canonical R2\n")
    _write_fastq(fastq, "S1_L001_R1.fastq.gz", "lane R1\n")
    _write_fastq(fastq, "S1_r2.fastq.gz", "lower R2\n")
    before = {p.name: p.read_bytes() for p in fastq.iterdir()}

    cfg = _settings(["S1"])
    assert mp.prepare_lane_split_fastqs(cfg) == 0
    assert {p.name: p.read_bytes() for p in fastq.iterdir()} == before
    assert mp.get_samples(cfg) == ["Test/S1_R1", "Test/S1_R2"]


def test_prepare_does_not_cross_sample_prefixes(fastq_cwd):
    fastq = fastq_cwd / "FASTQ" / "Test"
    for sample in ("S1", "S10"):
        _write_fastq(fastq, f"{sample}_L001_R1.fastq.gz", f"{sample} L001 R1\n")
        _write_fastq(fastq, f"{sample}_L002_R1.fastq.gz", f"{sample} L002 R1\n")
        _write_fastq(fastq, f"{sample}_R2.fastq.gz", f"{sample} R2\n")

    cfg = _settings(["S1", "S10"])
    assert mp.prepare_lane_split_fastqs(cfg) == 2
    assert _read_fastq(fastq / "S1_R1.fastq.gz") == "S1 L001 R1\nS1 L002 R1\n"
    assert _read_fastq(fastq / "S10_R1.fastq.gz") == "S10 L001 R1\nS10 L002 R1\n"
    assert mp.get_samples(cfg) == [
        "Test/S1_R1",
        "Test/S1_R2",
        "Test/S10_R1",
        "Test/S10_R2",
    ]


def test_prepare_orders_by_parsed_lane_and_chunk(fastq_cwd):
    fastq = fastq_cwd / "FASTQ" / "Test"
    _write_fastq(fastq, "S1_L010_R1.fastq.gz", "L010\n")
    _write_fastq(fastq, "S1_L002_R1_002.fastq.gz", "L002 chunk2\n")
    _write_fastq(fastq, "S1_L002_R1_001.fastq.gz", "L002 chunk1\n")
    _write_fastq(fastq, "S1_R2.fastq.gz", "R2\n")

    assert mp.prepare_lane_split_fastqs(_settings(["S1"])) == 1
    assert _read_fastq(fastq / "S1_R1.fastq.gz") == "L002 chunk1\nL002 chunk2\nL010\n"


def test_prepare_is_atomic_on_error(fastq_cwd, monkeypatch):
    fastq = fastq_cwd / "FASTQ" / "Test"
    _write_fastq(fastq, "S1_L001_R1.fastq.gz", "a\n")
    _write_fastq(fastq, "S1_L002_R1.fastq.gz", "b\n")
    _write_fastq(fastq, "S1_R2.fastq.gz", "c\n")

    calls = []
    real = shutil.copyfileobj

    def boom(src, dst, *args, **kwargs):
        calls.append(1)
        if len(calls) > 1:
            raise OSError("boom")
        return real(src, dst, *args, **kwargs)

    monkeypatch.setattr(mp.shutil, "copyfileobj", boom)
    assert mp.prepare_lane_split_fastqs(_settings(["S1"])) == 0
    assert not (fastq / "S1_R1.fastq.gz").exists()
    assert not (fastq / "S1_R1.fastq.gz.tmp").exists()


def test_six_raw_pairs(workdir):
    cfg, _fastq = workdir
    mp.prepare_lane_split_fastqs(cfg)
    samples = mp.get_samples(cfg)
    assert set(samples) == {f"Test/{s}_{r}" for s in SAMPLES for r in READS}
    assert not [x for x in samples if re.search(r"_L\d{1,3}", x)]


def test_generated_workflow_sample_lookup(workdir):
    cfg, _fastq = workdir
    conditions = mp.get_conditions(cfg)
    samples = mp.get_samples(cfg)
    combos = mw.get_combo(["QC"], cfg, conditions)
    subdir = "PreFlows"
    mp.create_skeleton(subdir, None)
    mw.make_pre("QC", cfg, samples, conditions, subdir, "INFO", "Pre", combinations=combos)

    generated = [
        p
        for p in glob.glob(os.path.join(subdir, "*_subsnake.smk"))
        if "def raw_fastq_input" in Path(p).read_text()
    ]
    assert generated
    text = Path(generated[0]).read_text()
    match = re.search(
        r"^def raw_fastq_input\(wildcards, read=None\):\n(?:.*\n)*?(?=^def )", text, re.M
    )
    assert match

    def _raw_fastq_input(samples, paired="paired"):
        namespace = {
            "os": os,
            "SAMPLES": samples,
            "paired": paired,
            "sample_from_wildcard": mp.sample_from_wildcard,
        }
        exec(match.group(0), namespace)
        return namespace["raw_fastq_input"]

    raw_fastq_input = _raw_fastq_input(
        ["Test/S1", "Test/S10", "Test/S2", "Test/S3", "Test/S8", "Test/S9"]
    )
    for sample in ("S1", "S10"):
        wildcards = types.SimpleNamespace(file=sample, rawfile=sample, read=None)
        assert raw_fastq_input(wildcards, "R1") == f"FASTQ/Test/{sample}_R1.fastq.gz"

    qualified = _raw_fastq_input(["CondA/S1", "CondB/S1"])
    for condition in ("CondA", "CondB"):
        for read in READS:
            wildcards = types.SimpleNamespace(
                file=f"{condition}/S1", rawfile=f"{condition}/S1", read=None
            )
            assert qualified(wildcards, read) == f"FASTQ/{condition}/S1_{read}.fastq.gz"
    assert qualified(
        types.SimpleNamespace(file="CondA/S1", rawfile="CondA/S1", read=None), "R1"
    ) != qualified(
        types.SimpleNamespace(file="CondB/S1", rawfile="CondB/S1", read=None), "R1"
    )

    for order in (["CondA/S1", "CondB/S1"], ["CondB/S1", "CondA/S1"]):
        ambiguous = _raw_fastq_input(order)
        with pytest.raises(ValueError) as excinfo:
            ambiguous(types.SimpleNamespace(file="S1", rawfile="S1", read=None), "R1")
        message = str(excinfo.value)
        assert "CondA/S1" in message and "CondB/S1" in message
        assert "condition qualified" in message

    single = _raw_fastq_input(["CondA/S1"])
    with pytest.raises(ValueError):
        single(
            types.SimpleNamespace(file="CondB/S1", rawfile="CondB/S1", read=None), "R1"
        )
    assert (
        single(types.SimpleNamespace(file="S1", rawfile="S1", read=None), "R1")
        == "FASTQ/CondA/S1_R1.fastq.gz"
    )

    repeated = _raw_fastq_input(["CondA/S1", "CondA/S1"])
    assert (
        repeated(types.SimpleNamespace(file="S1", rawfile="S1", read=None), "R1")
        == "FASTQ/CondA/S1_R1.fastq.gz"
    )

    singlecell = _raw_fastq_input(["Test/S1"], paired="singlecell")
    for read in READS:
        assert (
            singlecell(types.SimpleNamespace(file="S1", rawfile="S1", read=None), read)
            == f"FASTQ/Test/S1_{read}.fastq.gz"
        )


SELECTOR = re.compile(r"element -> return (.+)$", re.M)
SELECTOR_ELEMENTS = ("S1", "S10")


def _nf_post_selectors(cfg):
    conditions = mp.get_conditions(cfg)
    _pre, sub, _post = mw.get_processes(cfg)
    combos = mw.get_combo(sub, cfg, conditions)
    selectors = []
    for subwork in ("COUNTING", "CIRCS", "FUSIONS"):
        samples = mp.get_samples_postprocess(cfg, subwork)
        subdir = f"Post{subwork}"
        mp.create_skeleton(subdir, None)
        mw.nf_make_post(
            subwork, cfg, samples, conditions, subdir, "INFO", combinations=combos
        )
        for path in sorted(glob.glob(os.path.join(subdir, "*_subflow.nf"))):
            for match in SELECTOR.finditer(Path(path).read_text()):
                if match.group(1) not in selectors:
                    selectors.append(match.group(1))
    return selectors


def _resolve_selectors(engine, root, selectors):
    lines = ["params.gCOMBO = 'combo'", "workflow {"]
    for index, selector in enumerate(selectors):
        for element in SELECTOR_ELEMENTS:
            lines.append(
                f"    Channel.fromPath(['{element}'].collect{{ element -> return {selector} }})"
                f".view {{ it -> '{index}\\t{element}\\t' + it.name }}"
            )
    lines.append("}")
    workflow = root / "selectors.nf"
    workflow.write_text("\n".join(lines) + "\n")
    result = subprocess.run(
        [engine, "run", str(workflow), "-work-dir", str(root / "work")],
        capture_output=True,
        text=True,
        timeout=180,
        env={**os.environ, "JAVA_TOOL_OPTIONS": "-XX:-UseContainerSupport"},
    )
    assert result.returncode == 0, result.stdout + result.stderr
    resolved = {}
    for line in result.stdout.splitlines():
        match = re.match(r"^(\d+)\t(\S+)\t(\S+)$", line)
        if match:
            key = (int(match.group(1)), match.group(2))
            resolved.setdefault(key, set()).add(match.group(3))
    return resolved


def test_nextflow_selectors_do_not_cross_sample_prefixes(workdir, tmp_path):
    engine = shutil.which("nextflow")
    if not engine:
        pytest.skip("nextflow not on PATH")
    cfg, _fastq = workdir
    cfg["SETTINGS"]["Test"]["SEQUENCING"] = "single"
    cfg["CIRCS"] = {
        "TOOLS": {"ciri2": "CIRI2"},
        "Test": {
            "ciri2": {
                "OPTIONS": {
                    "CIRCANNO": "GENOMES/anno.gtf.gz",
                    "CIRCREF": "GENOMES/genome.fa.gz",
                }
            }
        },
    }
    cfg["FUSIONS"] = json.loads((DATA / "config_DTU_test.json").read_text())["FUSIONS"]

    mapped = tmp_path / "MAPPED" / "combo"
    trimmed = tmp_path / "TRIMMED_FASTQ" / "combo"
    dedup = tmp_path / "DEDUP_FASTQ" / "combo"
    for directory in (mapped, trimmed, dedup):
        directory.mkdir(parents=True)
    bam_variants = ("", "_unique", "_dedup", "_unique_dedup")
    for element in SELECTOR_ELEMENTS:
        for variant in bam_variants:
            (mapped / f"{element}_mapped_sorted{variant}.bam").touch()
        (mapped / f"{element}_mapped_sorted.sam.gz").touch()
        (mapped / f"{element}.Chimeric.out.junction").touch()
        (trimmed / f"{element}_trimmed.fastq.gz").touch()
        (dedup / f"{element}_dedup.fastq.gz").touch()

    selectors = _nf_post_selectors(cfg)
    resolved = _resolve_selectors(engine, tmp_path, selectors)

    def _hits(suffix):
        indexes = [
            index
            for index, selector in enumerate(selectors)
            if f'+element+"{suffix}"' in selector
        ]
        assert len(indexes) == 1, (suffix, indexes)
        return {
            element: resolved.get((indexes[0], element), set())
            for element in SELECTOR_ELEMENTS
        }

    assert _hits("_trimmed.fastq.gz") == {
        "S1": {"S1_trimmed.fastq.gz"},
        "S10": {"S10_trimmed.fastq.gz"},
    }
    assert _hits("_dedup.fastq.gz") == {
        "S1": {"S1_dedup.fastq.gz"},
        "S10": {"S10_dedup.fastq.gz"},
    }
    assert _hits("_mapped_sorted{,_*}.bam") == {
        element: {f"{element}_mapped_sorted{variant}.bam" for variant in bam_variants}
        for element in SELECTOR_ELEMENTS
    }
    assert _hits("_mapped_sorted.sam.gz") == {
        "S1": {"S1_mapped_sorted.sam.gz"},
        "S10": {"S10_mapped_sorted.sam.gz"},
    }
    assert _hits(".Chimeric.out.junction") == {
        "S1": {"S1.Chimeric.out.junction"},
        "S10": {"S10.Chimeric.out.junction"},
    }


def _resolve_fastq_patterns(engine, root, subflow, samples, paired):
    text = Path(subflow).read_text()
    match = re.search(
        r"^def fastq_patterns\(samples, mode\)\{\n(?:.*\n)*?\}\n", text, re.M
    )
    assert match, subflow
    script = (
        "nextflow.enable.dsl=2\n"
        f"params.gPAIRED = '{paired}'\n"
        + match.group(0)
        + "\nworkflow {\n"
        + f"    Channel.fromPath(fastq_patterns({samples!r}, 'paired')).view {{ it -> it.name }}\n"
        + "}\n"
    )
    workflow = root / "patterns.nf"
    workflow.write_text(script)
    result = subprocess.run(
        [engine, "run", str(workflow), "-work-dir", str(root / "work")],
        capture_output=True,
        text=True,
        timeout=180,
        env={**os.environ, "JAVA_TOOL_OPTIONS": "-XX:-UseContainerSupport"},
    )
    assert result.returncode == 0, result.stdout + result.stderr
    return {
        line.strip()
        for line in result.stdout.splitlines()
        if line.strip().endswith(".fastq.gz")
    }


def test_singlecell_nextflow_fastq_patterns_bind_canonical_mates(
    tmp_path, monkeypatch, repo_paths
):
    engine = shutil.which("nextflow")
    if not engine:
        pytest.skip("nextflow not on PATH")
    fastq = tmp_path / "FASTQ" / "Test"
    _touch(fastq, "S1_R1.fastq.gz", "S1_R2.fastq.gz")
    cfg = _qc_config()
    cfg["SETTINGS"]["Test"]["SAMPLES"] = ["S1"]
    cfg["SETTINGS"]["Test"]["SEQUENCING"] = "singlecell"
    configfile = tmp_path / "config_qc.json"
    configfile.write_text(json.dumps(cfg))
    monkeypatch.chdir(tmp_path)

    assert mp.get_samples(cfg) == ["Test/S1"]
    mr.run_nextflow(
        str(configfile), ".", True, 2, False, "INFO", save=True, optionalargs=[]
    )
    subflows = glob.glob(os.path.join("SubFlows", "*_subflow.nf"))
    assert subflows
    resolved = _resolve_fastq_patterns(
        engine, tmp_path, subflows[0], ["Test/S1"], "singlecell"
    )
    assert resolved == {"S1_R1.fastq.gz", "S1_R2.fastq.gz"}


def _settings(samples, sequencing="paired", condition=("Test",)):
    node = {
        "SAMPLES": list(samples),
        "SEQUENCING": sequencing,
        "GROUPS": [],
        "TYPES": "",
        "BATCHES": "",
        "INDEX": "",
        "PREFIX": "",
    }
    settings = {}
    current = settings
    for part in condition[:-1]:
        current = current.setdefault(part, {})
    current[condition[-1]] = node
    return {"SETTINGS": settings}


def _touch(directory, *names):
    directory.mkdir(parents=True, exist_ok=True)
    for name in names:
        (directory / name).write_bytes(b"")


@pytest.fixture
def fastq_cwd(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    yield tmp_path


def test_condition_parent_fallback_single_enumeration(fastq_cwd, monkeypatch):
    fastq = fastq_cwd / "FASTQ" / "Ecoli"
    _touch(
        fastq,
        "S1_R1.fastq.gz",
        "S1_R2.fastq.gz",
        "S2_R1.fastq.gz",
        "S2_R2.fastq.gz",
    )
    cfg = _settings(["S1", "S2"], condition=("Ecoli", "WT"))
    scans = []
    real_listdir = os.listdir
    monkeypatch.setattr(
        os, "listdir", lambda p=".": scans.append(str(p)) or real_listdir(p)
    )
    assert mp.get_samples_from_dir(["Ecoli", "WT"], cfg) == ["Ecoli/S1", "Ecoli/S2"]
    assert scans == [
        os.path.join("FASTQ", "Ecoli", "WT"),
        os.path.join("FASTQ", "Ecoli"),
    ]


def test_missing_sample_not_matched_by_longer_name(fastq_cwd):
    fastq = fastq_cwd / "FASTQ" / "Test"
    _touch(fastq, "S10_R1.fastq.gz", "S10_R2.fastq.gz")
    cfg = _settings(["S1", "S10"])
    assert mp.get_samples_from_dir(["Test"], cfg) == ["Test/S10"]
    assert mp.get_samples(cfg) == ["Test/S10_R1", "Test/S10_R2"]
    assert mp.check_samples(cfg) is False
    _touch(fastq, "S1_R1.fastq.gz", "S1_R2.fastq.gz")
    assert mp.check_samples(cfg) is True
    assert mp.get_samples(cfg) == [
        "Test/S1_R1",
        "Test/S1_R2",
        "Test/S10_R1",
        "Test/S10_R2",
    ]


def test_check_samples_without_configured_samples(fastq_cwd):
    _touch(fastq_cwd / "FASTQ" / "Test", "S1_R1.fastq.gz")
    assert mp.check_samples(_settings([])) is True


def test_sequencing_shapes(fastq_cwd):
    _touch(
        fastq_cwd / "FASTQ" / "Test",
        "S1_R1.fastq.gz",
        "S1_R2.fastq.gz",
    )
    paired = _settings(["S1"], sequencing="paired")
    assert mp.get_samples(paired) == ["Test/S1_R1", "Test/S1_R2"]
    assert mp.get_samples_from_dir(["Test"], paired) == ["Test/S1"]

    single = _settings(["S1"], sequencing="single")
    assert mp.get_samples(single) == ["Test/S1"]
    assert mp.get_samples_from_dir(["Test"], single) == ["Test/S1_R1", "Test/S1_R2"]

    singlecell = _settings(["S1"], sequencing="singlecell")
    assert mp.get_samples(singlecell) == ["Test/S1"]
    assert mp.get_samples_from_dir(["Test"], singlecell) == ["Test/S1"]


def test_merged_canonical_files_preferred(fastq_cwd):
    fastq = fastq_cwd / "FASTQ" / "Test"
    _touch(
        fastq,
        "S1_L001_R1_001.fastq.gz",
        "S1_L002_R1_001.fastq.gz",
        "S1_R1.fastq.gz",
        "S1_R2.fastq.gz",
    )
    cfg = _settings(["S1"])
    assert mp.get_samples(cfg) == ["Test/S1_R1", "Test/S1_R2"]
    assert mp.get_samples_from_dir(["Test"], cfg) == ["Test/S1"]


def test_discovery_does_not_rename_files(fastq_cwd):
    fastq = fastq_cwd / "FASTQ" / "Test"
    _touch(fastq, "S1_r1.fastq.gz", "S1_r2.fastq.gz")
    before = sorted(p.name for p in fastq.iterdir())
    cfg = _settings(["S1"])
    assert mp.get_samples_from_dir(["Test"], cfg) == ["Test/S1"]
    assert mp.get_samples(cfg) == ["Test/S1_r1", "Test/S1_r2"]
    assert sorted(p.name for p in fastq.iterdir()) == before


def test_discovery_cache_is_per_call(fastq_cwd):
    fastq = fastq_cwd / "FASTQ" / "Test"
    _touch(fastq, "S1_R1.fastq.gz", "S1_R2.fastq.gz")
    assert mp.get_samples(_settings(["S1"])) == ["Test/S1_R1", "Test/S1_R2"]
    _touch(fastq, "S2_R1.fastq.gz", "S2_R2.fastq.gz")
    cfg = _settings(["S1", "S2"])
    assert mp.get_samples(cfg) == [
        "Test/S1_R1",
        "Test/S1_R2",
        "Test/S2_R1",
        "Test/S2_R2",
    ]
    assert mp.check_samples(cfg) is True


def _measure_scans(tmp_path, monkeypatch, sample_count):
    root = tmp_path / f"n{sample_count}"
    names = []
    for index in range(sample_count):
        names.extend([f"S{index}_R1.fastq.gz", f"S{index}_R2.fastq.gz"])
    _touch(root / "FASTQ" / "Test", *names)
    cfg = _settings([f"S{index}" for index in range(sample_count)])
    scans = []
    matches = []
    real_listdir = os.listdir
    real_match = mp._matches_sample_fastq
    with monkeypatch.context() as patch:
        patch.chdir(root)
        patch.setattr(
            os, "listdir", lambda p=".": scans.append(str(p)) or real_listdir(p)
        )
        patch.setattr(
            mp,
            "_matches_sample_fastq",
            lambda f, s: matches.append((f, s)) or real_match(f, s),
        )
        samples = mp.get_samples(cfg)
    assert len(samples) == 2 * sample_count
    return scans, matches, len(names)


def test_directory_scans_do_not_grow_with_sample_count(tmp_path, monkeypatch):
    small_scans, small_matches, small_files = _measure_scans(tmp_path, monkeypatch, 2)
    large_scans, large_matches, large_files = _measure_scans(
        tmp_path, monkeypatch, 40
    )
    assert small_scans == large_scans == [os.path.join("FASTQ", "Test")]
    assert len(small_matches) == small_files
    assert len(large_matches) <= 2 * large_files
    assert len(large_matches) < large_files * 40


@pytest.fixture
def repo_paths(monkeypatch):
    monkeypatch.setattr(mw, "workflowpath", str(REPO / "workflows"))
    monkeypatch.setattr(mw, "envpath", str(REPO / "envs") + os.sep)


def _qc_config():
    cfg = _settings(["S1", "S10"])
    cfg["SETTINGS"]["Test"].update(
        {"REFERENCE": "", "DECOY": "", "ANNOTATION": {"GTF": "", "GFF": ""}, "IP": ""}
    )
    cfg.update(
        {
            "WORKFLOWS": "QC",
            "BINS": "",
            "MAXTHREADS": "2",
            "VERSION": MONSDA.__version__,
            "QC": {
                "TOOLS": {"fastqc": "fastqc"},
                "Test": {"fastqc": {"OPTIONS": {"QC": "", "MULTI": ""}}},
            },
        }
    )
    return cfg


@pytest.mark.parametrize("engine", ["snakemake", "nextflow"])
def test_runner_keeps_qc_step_for_fastq_only_config(
    tmp_path, monkeypatch, repo_paths, engine
):
    fastq = tmp_path / "FASTQ" / "Test"
    _touch(fastq, "S1_R1.fastq.gz", "S1_R2.fastq.gz", "S10_R1.fastq.gz", "S10_R2.fastq.gz")
    configfile = tmp_path / "config_qc.json"
    configfile.write_text(json.dumps(_qc_config()))
    monkeypatch.chdir(tmp_path)

    if engine == "snakemake":
        mr.run_snakemake(
            str(configfile), ".", True, 2, False, "INFO", save=True, optionalargs=[]
        )
        expected = tmp_path / "SubSnakes" / "Test_PreQC_fastqc_subsnake.smk"
    else:
        mr.run_nextflow(
            str(configfile), ".", True, 2, False, "INFO", save=True, optionalargs=[]
        )
        expected = tmp_path / "SubFlows" / "Test_pre_QC_fastqc_subflow.nf"

    commands = (tmp_path / "JOBS" / "MONSDA.commands").read_text().splitlines()
    commands = [line for line in commands if line.strip()]
    assert len(commands) == 1, commands
    assert expected.name in commands[0], commands[0]
    assert expected.is_file()
    if engine == "snakemake":
        confo = tmp_path / "SubSnakes" / "Test_PreQC_fastqc_subconfig.json"
        assert json.loads(confo.read_text())["SAMPLES"] == ["S1", "S10"], confo
    else:
        assert "Test/S1" in commands[0] and "Test/S10" in commands[0], commands[0]


BASECALL_RAW = {"dorado": ".pod5", "guppy": ".fast5"}
BASECALL_SEQ = {"S1": "ACGTACGT", "S10": "TTTTGGGG"}
BASECALL_FOREIGN = {"S1": "GGGGCCCC", "S10": "CCCCGGGG"}


@pytest.fixture
def rawdir(tmp_path, monkeypatch):
    raw = tmp_path / "RAW" / "Test"
    raw.mkdir(parents=True)
    for suffix in BASECALL_RAW.values():
        for sample in BASECALL_SEQ:
            (raw / f"{sample}{suffix}").write_text(BASECALL_SEQ[sample] + "\n")
    monkeypatch.chdir(tmp_path)
    yield tmp_path


def _basecall_config(tool):
    cfg = _settings(["S1", "S10"], sequencing="single")
    cfg["SETTINGS"]["Test"].update(
        {"REFERENCE": "", "DECOY": "", "ANNOTATION": {"GTF": "", "GFF": ""}, "IP": ""}
    )
    cfg.update(
        {
            "WORKFLOWS": "BASECALL",
            "BINS": "",
            "MAXTHREADS": "2",
            "VERSION": MONSDA.__version__,
            "BASECALL": {
                "TOOLS": {tool: tool},
                "Test": {tool: {"OPTIONS": {"BASECALL": "", "MODEL": "model"}}},
            },
        }
    )
    return cfg


def test_raw_signal_map_exact_files_and_foreign_format(rawdir):
    cfg = _basecall_config("guppy")
    assert mp.basecall_signal_map(cfg) == {
        "Test/S1": ["RAW/Test/S1.fast5"],
        "Test/S10": ["RAW/Test/S10.fast5"],
    }
    assert mp.basecall_signal_map(cfg, "dorado") == {
        "Test/S1": ["RAW/Test/S1.pod5"],
        "Test/S10": ["RAW/Test/S10.pod5"],
    }
    assert mp.basecall_samples(cfg) == ["Test/S1", "Test/S10"]


def test_raw_signal_map_recursive_directory(rawdir):
    raw = rawdir / "RAW" / "Test"
    (raw / "S1.fast5").unlink()
    for sub, seq in (("a", "AAAA"), ("b", "CCCC")):
        (raw / "S1" / sub).mkdir(parents=True)
        (raw / "S1" / sub / "run.fast5").write_text(seq + "\n")
    (raw / "S1" / "notes.txt").write_text("ignored\n")
    assert mp.basecall_signal_map(_basecall_config("guppy")) == {
        "Test/S1": ["RAW/Test/S1/a/run.fast5", "RAW/Test/S1/b/run.fast5"],
        "Test/S10": ["RAW/Test/S10.fast5"],
    }


def test_raw_signal_map_explicit_list(rawdir):
    raw = rawdir / "RAW" / "Test"
    (raw / "S1.fast5").unlink()
    (raw / "loose").mkdir()
    (raw / "loose" / "one.fast5").write_text("AAAA\n")
    (raw / "loose" / "two.fast5").write_text("CCCC\n")
    (raw / "S1.fast5.list").write_text("loose/one.fast5\nloose/two.fast5\n")
    assert mp.basecall_signal_map(_basecall_config("guppy")) == {
        "Test/S1": ["RAW/Test/loose/one.fast5", "RAW/Test/loose/two.fast5"],
        "Test/S10": ["RAW/Test/S10.fast5"],
    }


def test_raw_signal_map_does_not_prefix_match(rawdir):
    raw = rawdir / "RAW" / "Test"
    (raw / "S1.fast5").unlink()
    (raw / "S1_extra.fast5").write_text("AAAA\n")
    with pytest.raises(ValueError) as excinfo:
        mp.basecall_signal_map(_basecall_config("guppy"))
    assert "no raw signal source for sample Test/S1" in str(excinfo.value)


def test_raw_signal_map_accepts_symlinked_signal(rawdir):
    raw = rawdir / "RAW" / "Test"
    (raw / "S1.fast5").unlink()
    (raw / "real").mkdir()
    (raw / "real" / "one.fast5").write_text("AAAA\n")
    (raw / "S1.fast5").symlink_to(raw / "real" / "one.fast5")
    assert mp.basecall_signal_map(_basecall_config("guppy")) == {
        "Test/S1": ["RAW/Test/S1.fast5"],
        "Test/S10": ["RAW/Test/S10.fast5"],
    }


def test_raw_signal_map_rejects_symlinked_duplicate(rawdir):
    raw = rawdir / "RAW" / "Test"
    (raw / "S1.fast5").unlink()
    (raw / "alias.fast5").symlink_to(raw / "S10.fast5")
    (raw / "S1.fast5.list").write_text("alias.fast5\n")
    with pytest.raises(ValueError) as excinfo:
        mp.basecall_signal_map(_basecall_config("guppy"))
    assert "is claimed by samples" in str(excinfo.value)


def test_raw_signal_map_requires_selected_caller(rawdir):
    cfg = _basecall_config("guppy")
    cfg["BASECALL"]["TOOLS"] = {"guppy": "guppy", "dorado": "dorado"}
    assert mp.basecall_caller(cfg) is None
    with pytest.raises(ValueError) as excinfo:
        mp.basecall_signal_map(cfg)
    assert "no basecaller selected" in str(excinfo.value)
    assert mp.basecall_samples(cfg) == ["Test/S1", "Test/S10"]


def _break_source(raw, mutation, suffix):
    if mutation == "missing":
        (raw / f"S1{suffix}").unlink()
    elif mutation == "empty":
        (raw / f"S1{suffix}").write_bytes(b"")
    elif mutation == "ambiguous":
        (raw / "S1").mkdir()
        (raw / "S1" / f"run{suffix}").write_text("AAAA\n")
    elif mutation == "emptydir":
        (raw / f"S1{suffix}").unlink()
        (raw / "S1").mkdir()
        (raw / "S1" / "notes.txt").write_text("ignored\n")
    elif mutation == "emptylist":
        (raw / f"S1{suffix}").unlink()
        (raw / f"S1{suffix}.list").write_text("\n")
    elif mutation == "overlap":
        (raw / f"S1{suffix}").unlink()
        (raw / f"S1{suffix}.list").write_text(f"S10{suffix}\n")
    elif mutation == "duplicate":
        (raw / f"S1{suffix}").unlink()
        (raw / f"S1{suffix}.list").write_text(f"S10{suffix}\nS10{suffix}\n")
    elif mutation == "outside":
        (raw / f"S1{suffix}").unlink()
        (raw / f"S1{suffix}.list").write_text(f"../Other/run{suffix}\n")
    elif mutation == "absolute":
        (raw / f"S1{suffix}").unlink()
        (raw / f"S1{suffix}.list").write_text(f"/tmp/run{suffix}\n")
    elif mutation == "wrongtype":
        (raw / f"S1{suffix}").unlink()
        other = [x for x in BASECALL_RAW.values() if x != suffix][0]
        (raw / f"S1{suffix}.list").write_text(f"S10{other}\n")
    elif mutation == "notregular":
        (raw / f"S1{suffix}").unlink()
        (raw / "loose").mkdir()
        (raw / f"S1{suffix}.list").write_text("loose\n")
    else:
        raise AssertionError(mutation)


@pytest.mark.parametrize(
    "mutation,message",
    [
        ("missing", "no raw signal source for sample Test/S1"),
        ("empty", "is empty"),
        ("ambiguous", "ambiguous raw signal source for sample Test/S1"),
        ("emptydir", "contains no .fast5 files"),
        ("emptylist", "names no raw signal file"),
        ("overlap", "is claimed by samples"),
        ("duplicate", "duplicate entry"),
        ("outside", "points outside"),
        ("absolute", "is absolute"),
        ("wrongtype", "does not end in .fast5"),
        ("notregular", "is not a regular file"),
    ],
)
def test_raw_signal_map_rejects_invalid_sources(rawdir, mutation, message):
    _break_source(rawdir / "RAW" / "Test", mutation, ".fast5")
    with pytest.raises(ValueError) as excinfo:
        mp.basecall_signal_map(_basecall_config("guppy"))
    assert message in str(excinfo.value), str(excinfo.value)


STUB_GUPPY = """#!/bin/bash
printf '%s\\n' "$*" >> "$STUB_LOG"
indir="."; list=""; out=""
while [ $# -gt 0 ]; do
  case "$1" in
    -i) indir="$2"; shift 2;;
    --input_file_list) list="$2"; shift 2;;
    -s) out="$2"; shift 2;;
    *) shift;;
  esac
done
mkdir -p "$out/pass"
: > "$out/pass/records.fastq"
i=0
while IFS= read -r n; do
  [ -z "$n" ] && continue
  [ -f "$indir/$n" ] || { echo "missing $indir/$n" >&2; exit 1; }
  i=$((i+1))
  printf '@read%d\\n%s\\n+\\n+\\n' "$i" "$(cat "$indir/$n")" >> "$out/pass/records.fastq"
done < "$list"
gzip -c "$out/pass/records.fastq" > "$out/pass/fastq_runid_1.fastq.gz"
rm -f "$out/pass/records.fastq"
echo summary > "$out/sequencing_summary.txt"
echo telemetry > "$out/sequencing_telemetry.js"
echo log > "$out/guppy.log"
"""

STUB_DORADO = """#!/bin/bash
printf '%s\\n' "$*" >> "$STUB_LOG"
cmd="$1"; shift
if [ "$cmd" = "download" ]; then
  model=""; dir="."
  while [ $# -gt 0 ]; do
    case "$1" in
      --model) model="$2"; shift 2;;
      --directory) dir="$2"; shift 2;;
      *) shift;;
    esac
  done
  mkdir -p "$dir/$model"
  echo model > "$dir/$model/model.stamp"
  exit 0
fi
[ "$cmd" = "basecaller" ] || { echo "unexpected command $cmd" >&2; exit 1; }
[ $# -eq 2 ] || { echo "expected model and one data directory, got $#: $*" >&2; exit 1; }
model="$1"; data="$2"
[ -d "$model" ] || { echo "model $model is not a directory" >&2; exit 1; }
[ -d "$data" ] || { echo "data $data is not a directory" >&2; exit 1; }
reads=$(mktemp)
trap 'rm -f "$reads"' EXIT
for f in "$data"/*; do cat "$f" >> "$reads"; done
printf '@HD\\tVN:1.6\\n'
i=0
while IFS= read -r line; do
  i=$((i+1))
  printf 'read%d\\t4\\t*\\t0\\t0\\t*\\t*\\t0\\t0\\t%s\\t*\\n' "$i" "$line"
done < "$reads"
"""


@pytest.fixture
def basecall_stubs(tmp_path, monkeypatch):
    stubs = tmp_path / "stubs"
    stubs.mkdir()
    for name, body in (("guppy", STUB_GUPPY), ("dorado", STUB_DORADO)):
        path = stubs / name
        path.write_text(body)
        path.chmod(0o755)
    monkeypatch.setenv("PATH", str(stubs) + os.pathsep + os.environ["PATH"])
    monkeypatch.setenv("STUB_LOG", str(tmp_path / "stub.log"))
    yield tmp_path / "stub.log"


def _require_engine(engine, tool):
    if not shutil.which("nextflow" if engine == "nf" else "snakemake"):
        pytest.skip(f"{engine} not on PATH")
    if tool == "dorado":
        for binary in ("samtools", "pigz"):
            if not shutil.which(binary):
                pytest.skip(f"{binary} not on PATH")


def _run_basecall(engine, tool, root):
    cfg = _basecall_config(tool)
    conditions = mp.get_conditions(cfg)
    samples = mp.basecall_samples(cfg)
    subdir = f"Pre{engine}{tool}"
    mp.create_skeleton(subdir, None)
    if engine == "nf":
        jobs = mw.nf_make_pre("BASECALL", cfg, samples, conditions, subdir, "INFO")
        assert len(jobs) == 1, jobs
        nfo, _confo, tp, params = jobs[0]
        manifest = Path(params["BASECALL_MANIFEST"])
        assert manifest.is_file(), manifest
        assert manifest.name.endswith(".signals.json"), manifest
        assert manifest.parent == Path(nfo).parent, manifest
        assert json.loads(manifest.read_text()) == mp.basecall_signal_map(cfg, tool)
        pars = " ".join(f"--{key} {value}" for key, value in params.items())
        result = subprocess.run(
            f"nextflow -log {root}/run.log run {nfo} -w {root}/NextFlowWork {pars} {tp}",
            shell=True,
            capture_output=True,
            text=True,
            timeout=600,
            env={**os.environ, "JAVA_TOOL_OPTIONS": "-XX:-UseContainerSupport"},
        )
    else:
        jobs = mw.make_pre(
            "BASECALL", cfg, samples, conditions, subdir, "INFO", "Pre"
        )
        assert len(jobs) == 1, jobs
        smko, confo = jobs[0]
        result = subprocess.run(
            [
                shutil.which("snakemake"),
                "-j",
                "2",
                "-s",
                smko,
                "--configfile",
                confo,
                "--directory",
                ".",
                "--nolock",
            ],
            capture_output=True,
            text=True,
            timeout=600,
        )
    assert result.returncode == 0, result.stdout + result.stderr


def _fastq_sequences(root, sample):
    path = root / "FASTQ" / "Test" / f"{sample}.fastq.gz"
    assert path.is_file(), path
    with gzip.open(path, "rt") as fh:
        lines = fh.read().splitlines()
    assert len(lines) % 4 == 0, lines
    return [lines[index + 1] for index in range(0, len(lines), 4)]


def _layout_single(raw, suffix):
    return {"S1": [BASECALL_SEQ["S1"]], "S10": [BASECALL_SEQ["S10"]]}


def _layout_nested(raw, suffix):
    (raw / f"S1{suffix}").unlink()
    for sub, seq in (("a", "AAAA"), ("b", "CCCC")):
        (raw / "S1" / sub).mkdir(parents=True)
        (raw / "S1" / sub / f"run{suffix}").write_text(seq + "\n")
    return {"S1": ["AAAA", "CCCC"], "S10": [BASECALL_SEQ["S10"]]}


def _layout_list(raw, suffix):
    (raw / f"S1{suffix}").unlink()
    (raw / "loose").mkdir()
    (raw / "loose" / f"one{suffix}").write_text("AAAA\n")
    (raw / "loose" / f"two{suffix}").write_text("CCCC\n")
    (raw / f"S1{suffix}.list").write_text(f"loose/one{suffix}\nloose/two{suffix}\n")
    return {"S1": ["AAAA", "CCCC"], "S10": [BASECALL_SEQ["S10"]]}


BASECALL_LAYOUTS = {
    "single": _layout_single,
    "nested": _layout_nested,
    "list": _layout_list,
}


@pytest.mark.parametrize("engine", ["nf", "smk"])
@pytest.mark.parametrize("tool", sorted(BASECALL_RAW))
@pytest.mark.parametrize("layout", sorted(BASECALL_LAYOUTS))
def test_basecall_workflow_per_sample_output(
    rawdir, repo_paths, basecall_stubs, engine, tool, layout
):
    _require_engine(engine, tool)
    raw = rawdir / "RAW" / "Test"
    other = [x for x in BASECALL_RAW.values() if x != BASECALL_RAW[tool]][0]
    for sample, foreign in BASECALL_FOREIGN.items():
        (raw / f"{sample}{other}").write_text(foreign + "\n")
    expected = BASECALL_LAYOUTS[layout](raw, BASECALL_RAW[tool])
    before = {p: p.read_bytes() for p in raw.rglob("*") if p.is_file()}

    _run_basecall(engine, tool, rawdir)

    for sample, sequences in expected.items():
        assert sorted(_fastq_sequences(rawdir, sample)) == sorted(sequences)
    assert {p: p.read_bytes() for p in raw.rglob("*") if p.is_file()} == before
    for sample in BASECALL_SEQ:
        if tool == "guppy":
            assert (rawdir / "FASTQ" / "Test" / f"{sample}_summary.txt").is_file()
            assert (rawdir / "FASTQ" / "Test" / f"{sample}_telemetry.js").is_file()
        if engine == "nf":
            logpath = rawdir / "LOGS" / "Test" / "BASECALL" / tool / f"{sample}.log"
        else:
            logpath = (
                rawdir / "LOGS" / "Test" / sample / "BASECALL" / tool / f"{tool}.log"
            )
        assert logpath.is_file(), logpath


@pytest.mark.parametrize("engine", ["nf", "smk"])
@pytest.mark.parametrize("tool", sorted(BASECALL_RAW))
@pytest.mark.parametrize("mutation", ["missing", "ambiguous", "overlap"])
def test_basecall_workflow_rejects_invalid_sources(
    rawdir, repo_paths, engine, tool, mutation
):
    _break_source(rawdir / "RAW" / "Test", mutation, BASECALL_RAW[tool])
    with pytest.raises(SystemExit):
        mp.basecall_samples(_basecall_config(tool))


POST_STEPS = [
    "COUNTING",
    "TRACKS",
    "PEAKS",
    "DE",
    "DEU",
    "DAS",
    "DTU",
    "CIRCS",
    "FUSIONS",
]


def _full_config():
    cfg = json.loads((DATA / "config_Test.json").read_text())
    cfg["SETTINGS"]["Test"]["SAMPLES"] = ["S1", "S10"]
    cfg["SETTINGS"]["Test"]["GROUPS"] = ["T1", "T2"]
    cfg["WORKFLOWS"] = "QC,TRIMMING,DEDUP,MAPPING," + ",".join(POST_STEPS)
    cfg["FUSIONS"] = json.loads((DATA / "config_DTU_test.json").read_text())["FUSIONS"]
    cfg["CIRCS"] = {
        "TOOLS": {"ciri2": "CIRI2"},
        "Test": {
            "ciri2": {
                "OPTIONS": {
                    "CIRCANNO": "GENOMES/anno.gtf.gz",
                    "CIRCREF": "GENOMES/genome.fa.gz",
                }
            }
        },
    }
    return cfg


def _assert_unique(confos, expected):
    assert confos
    for confo in confos:
        ids = json.loads(Path(confo).read_text())["SAMPLES"]
        assert len(ids) == len(set(ids)), (confo, ids)
        assert sorted(ids) == sorted(expected), (confo, ids)


@pytest.mark.parametrize("engine", ["smk", "nf"])
def test_generated_configs_have_unique_sample_ids(fastq_cwd, repo_paths, engine):
    fastq = fastq_cwd / "FASTQ" / "Test"
    _touch(
        fastq,
        "S1_R1.fastq.gz",
        "S1_R2.fastq.gz",
        "S10_R1.fastq.gz",
        "S10_R2.fastq.gz",
    )

    configfile = fastq_cwd / "config_qc.json"
    configfile.write_text(json.dumps(_qc_config()))
    if engine == "smk":
        mr.run_snakemake(
            str(configfile), ".", True, 2, False, "INFO", save=True, optionalargs=[]
        )
        _assert_unique(glob.glob("SubSnakes/*_subconfig.json"), ["S1", "S10"])
    else:
        mr.run_nextflow(
            str(configfile), ".", True, 2, False, "INFO", save=True, optionalargs=[]
        )
        _assert_unique(glob.glob("SubFlows/*_subconfig.json"), ["S1", "S10"])

    cfg = _full_config()
    conditions = mp.get_conditions(cfg)
    _pre, sub, _post = mw.get_processes(cfg)
    combos = mw.get_combo(sub, cfg, conditions)
    samples = mp.get_samples(cfg)
    assert sorted(samples) == sorted(
        [
            "Test/S1_R1",
            "Test/S1_R2",
            "Test/S10_R1",
            "Test/S10_R2",
        ]
    )

    subdir = f"Sub{engine}"
    mp.create_skeleton(subdir, None)
    if engine == "smk":
        jobs = mw.make_sub(
            sub, cfg, samples, conditions, subdir, "INFO", combinations=combos
        )
    else:
        jobs = mw.nf_make_sub(
            sub, cfg, samples, conditions, subdir, "INFO", combinations=combos
        )
    _assert_unique([job[1] for job in jobs], ["S1", "S10"])

    for step in POST_STEPS:
        samples = mp.get_samples_postprocess(cfg, step)
        assert sorted(samples) == ["Test/S1", "Test/S10"], (step, samples)
        subdir = f"Post{engine}{step}"
        mp.create_skeleton(subdir, None)
        if engine == "smk":
            jobs = mw.make_post(
                step, cfg, samples, conditions, subdir, "INFO", combinations=combos
            )
        else:
            jobs = mw.nf_make_post(
                step, cfg, samples, conditions, subdir, "INFO", combinations=combos
            )
        _assert_unique([job[1] for job in jobs], ["Test/S1", "Test/S10"])
