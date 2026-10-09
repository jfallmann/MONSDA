"""Strandedness CLI and executable workflow gates; no writes to shared fixtures."""
import importlib.util
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
CHECKER = ROOT / 'scripts/Analysis/check_strandedness.py'
spec = importlib.util.spec_from_file_location('check_strandedness', CHECKER)
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)


@pytest.fixture(autouse=True)
def sandbox(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)


def infer_text(paired=True, forward=.95, reverse=.03, failed=.02):
    patterns = (checker.PAIRED_FORWARD, checker.PAIRED_REVERSE) if paired else (
        checker.SINGLE_FORWARD, checker.SINGLE_REVERSE)
    return (
        f'This is {"Pair" if paired else "Single"}End Data\n'
        f'Fraction of reads failed to determine: {failed:.4f}\n'
        f'Fraction of reads explained by "{patterns[0]}": {forward:.4f}\n'
        f'Fraction of reads explained by "{patterns[1]}": {reverse:.4f}\n'
    )


def run_checker(tmp_path, text, expected='forward', flat=False):
    qc = tmp_path / 'QC sample'
    folder = qc if flat else qc / 'rseqc/infer_experiment'
    folder.mkdir(parents=True, exist_ok=True)
    (folder / 'sample.infer_experiment.txt').write_text(text)
    report = qc / 'strandedness_check.json'
    result = subprocess.run(
        [sys.executable, str(CHECKER), '--qc-dir', str(qc), '--expected', expected,
         '--sample', 'sample name', '--report', str(report)],
        capture_output=True, text=True, timeout=20)
    return result, json.loads(report.read_text()) if report.exists() else None


@pytest.mark.parametrize('paired', [True, False])
@pytest.mark.parametrize('flat', [True, False])
@pytest.mark.parametrize('inferred,forward,reverse', [
    ('forward', .95, .03), ('reverse', .03, .95), ('unstranded', .49, .49)])
@pytest.mark.parametrize('expected', ['forward', 'reverse', 'unstranded', 'fr', 'rf'])
def test_verdict(tmp_path, paired, flat, inferred, forward, reverse, expected):
    result, report = run_checker(tmp_path, infer_text(paired, forward, reverse), expected, flat)
    match = checker.EXPECTED_ALIASES[expected] == inferred
    assert result.returncode == (0 if match else 2), result.stderr
    assert report['status'] == ('match' if match else 'mismatch')
    assert report['inferred'] == inferred
    assert report['paired'] is paired
    if not match:
        assert 'ERROR: strandedness mismatch for sample sample name' in result.stderr
        assert 'SEQUENCING' in result.stderr


@pytest.mark.parametrize('forward,reverse,failed,status', [
    (.2, .2, .6, 'inconclusive'), (.49, .21, .3, 'ambiguous'),
    (0, 0, 0, 'inconclusive'), (0, 0, 1, 'inconclusive'),
    (.3333, .3333, .3333, 'mismatch')])
def test_uncertainty_and_rounding(tmp_path, forward, reverse, failed, status):
    result, report = run_checker(tmp_path, infer_text(True, forward, reverse, failed))
    assert report['status'] == status
    assert result.returncode == (2 if status == 'mismatch' else 0)
    if status != 'mismatch':
        assert 'WARNING: strandedness NOT validated' in result.stderr


@pytest.mark.parametrize('forward,reverse,inferred', [
    (.8, .2, 'forward'), (.2, .8, 'reverse'), (.4, .6, 'unstranded'),
    (.6, .4, 'unstranded')])
def test_threshold_boundaries(forward, reverse, inferred):
    assert checker.infer_strandedness(0, forward, reverse) == inferred


@pytest.mark.parametrize('text', [
    '', infer_text() + 'Fraction of reads failed to determine: 0.0200\n',
    infer_text().replace('0.9500', 'NaN'), infer_text().replace('0.9500', 'inf'),
    infer_text().replace('0.9500', '-0.1'), infer_text().replace('0.9500', '1.1'),
    infer_text().replace('0.9500', 'bad'), infer_text().replace('0.9500', '0.1'),
    infer_text().replace('This is PairEnd Data\n', ''),
    infer_text().replace('Fraction of reads failed to determine: 0.0200\n', ''),
    infer_text().replace('1++,1--,2+-,2-+', 'unknown'),
    infer_text() + 'Fraction of reads explained by "++,--": 0.1\n',
    infer_text(False) + 'Fraction of reads explained by "1++,1--,2+-,2-+": 0.1\n',
])
def test_bad_report_fails(tmp_path, text):
    result, report = run_checker(tmp_path, text)
    assert result.returncode == 2
    assert 'ERROR:' in result.stderr
    assert report is None


def test_pattern_whitespace(tmp_path):
    result, report = run_checker(tmp_path, infer_text().replace(',', ', '))
    assert result.returncode == 0
    assert report['status'] == 'match'


def test_missing_and_multiple_reports(tmp_path):
    with pytest.raises(checker.StrandednessError, match='no .*found'):
        checker.find_infer_file(tmp_path)
    for name in ('a', 'b'):
        (tmp_path / f'{name}.infer_experiment.txt').write_text(infer_text())
    with pytest.raises(checker.StrandednessError, match='multiple'):
        checker.find_infer_file(tmp_path)


@pytest.fixture
def fake_rustqc(tmp_path, monkeypatch):
    """Only RustQC is stubbed; engines and the strand checker execute for real."""
    bin_dir = tmp_path / 'bin'
    bin_dir.mkdir()
    stub = bin_dir / 'rustqc'
    stub.write_text('#!/usr/bin/env python\n' +
        'import os, pathlib, sys\n'
        'args = sys.argv\n'
        'out = pathlib.Path(args[args.index("-o") + 1])\n'
        'out.mkdir(parents=True, exist_ok=True)\n'
        '(out / "rustqc_summary.json").write_text("{}")\n'
        '(out / "sample.infer_experiment.txt").write_text(os.environ["INFER_TEXT"])\n')
    stub.chmod(0o755)
    monkeypatch.setenv('PATH', str(bin_dir) + os.pathsep + os.environ.get('PATH', ''))
    shutil.copy2(ROOT / 'envs/rustqc.yaml', tmp_path / 'rustqc.yaml')
    bins = tmp_path / 'scripts with spaces'
    (bins / 'Analysis').mkdir(parents=True)
    shutil.copy2(CHECKER, bins / 'Analysis/check_strandedness.py')
    (tmp_path / 'annotation with spaces.gtf').write_text('fixture annotation\n')
    return bins


@pytest.mark.parametrize('unique', [False, True])
@pytest.mark.parametrize('mismatch', [False, True])
def test_snakemake_gate(tmp_path, monkeypatch, fake_rustqc, unique, mismatch):
    engine = shutil.which('snakemake')
    if not engine:
        pytest.skip('snakemake not on PATH')
    monkeypatch.setenv('INFER_TEXT', infer_text())
    bam = 'MAPPED/combo/sample_mapped_sorted' + ('_unique' if unique else '')
    (tmp_path / bam).parent.mkdir(parents=True)
    (tmp_path / (bam + '.bam')).touch()
    (tmp_path / (bam + '.bam.bai')).touch()
    out = bam.replace('MAPPED/', 'QC/')
    prefix = (
        "def env_bin_from_config(*args): return 'rustqc', 'rustqc'\n"
        "def tool_params(*args): return {'OPTIONS': {}}\n"
        "SAMPLES = ['sample']\nMAXTHREAD = 1\npaired = 'paired'\n"
        f"stranded = {'rf' if mismatch else 'fr'!r}\n"
        f"BINS = {str(fake_rustqc)!r}\nANNOTATION = 'annotation with spaces.gtf'\n"
    )
    downstream = (f'\nrule downstream:\n    input: "{out}/rustqc_summary.json"\n'
                  '    output: "downstream.done"\n    shell: "touch {output}"\n')
    snake = tmp_path / 'Snakefile'
    snake.write_text(prefix + (ROOT / 'workflows/rustqc.smk').read_text() + downstream)
    result = subprocess.run([engine, '-s', str(snake), '-j', '1', 'downstream.done'],
                            capture_output=True, text=True, timeout=90)
    assert (result.returncode != 0) is mismatch, result.stdout + result.stderr
    assert (tmp_path / 'downstream.done').exists() is not mismatch
    if mismatch:
        logs = '\n'.join(p.read_text() for p in (tmp_path / 'LOGS').rglob('*.log'))
        assert 'strandedness mismatch' in logs
    else:
        assert json.loads((tmp_path / out / 'strandedness_check.json').read_text())['status'] == 'match'


@pytest.mark.parametrize('mismatch', [False, True])
def test_nextflow_gate(tmp_path, monkeypatch, fake_rustqc, mismatch):
    engine = shutil.which('nextflow')
    if not engine:
        pytest.skip('nextflow not on PATH')
    from MONSDA.Workflows import nf_paramize
    monkeypatch.setenv('INFER_TEXT', infer_text())
    (tmp_path / 'sample.bam').touch()
    values = dict(POSTQCENV='rustqc', POSTQCBIN='rustqc', rustqc_params_QC='',
                  MAPPINGANNO='annotation with spaces.gtf')
    prefix = (
        'def get_always(key) {\n    return params[key]\n}\n'
        f'BINS = {json.dumps(str(fake_rustqc))}\n'
        f'STRANDED = {json.dumps("rf" if mismatch else "fr")}\n'
        "PAIRED = 'paired'\nTHREADS = 1\nCOMBO = 'combo'\nCONDITION = 'condition'\n"
    )
    suffix = '''
process downstream {
    input:
    path qc
    output:
    path 'downstream.done'
    script:
    "touch downstream.done"
}
workflow {
    rustqc_mapped(Channel.fromPath('sample.bam'))
    downstream(rustqc_mapped.out.rustqc_json)
}
'''
    workflow = tmp_path / 'test.nf'
    workflow.write_text(nf_paramize(prefix + (ROOT / 'workflows/rustqc.nf').read_text() + suffix))
    params = tmp_path / 'params.json'
    params.write_text(json.dumps(values))
    result = subprocess.run([engine, 'run', str(workflow), '-params-file', str(params),
                             '-work-dir', str(tmp_path / 'work')],
                            capture_output=True, text=True, timeout=120)
    assert (result.returncode != 0) is mismatch, result.stdout + result.stderr
    markers = list((tmp_path / 'work').rglob('downstream.done'))
    assert bool(markers) is not mismatch
    reports = list((tmp_path / 'work').rglob('strandedness_check.json'))
    assert reports, result.stdout + result.stderr
    assert json.loads(reports[0].read_text())['status'] == ('mismatch' if mismatch else 'match')
    if mismatch:
        assert 'strandedness mismatch' in result.stdout + result.stderr
        assert len(reports) == 1, 'Deterministic strand mismatch must not be retried'


@pytest.mark.parametrize('stream', ['stdout', 'stderr', 'silent'])
def test_runjob_failure_is_nonzero(stream):
    from MONSDA.RunMONSDA import runjob
    code = 'import sys; '
    if stream != 'silent':
        code += f'print("ERROR: strandedness mismatch", file=sys.{stream}); '
    code += 'sys.exit(2)'
    with pytest.raises(SystemExit) as exc:
        runjob(shlex.join([sys.executable, '-c', code]))
    assert exc.value.code not in (None, 0)


def test_runjob_success():
    from MONSDA.RunMONSDA import runjob
    assert runjob(shlex.join([sys.executable, '-c', 'print("completed")'])) == 0
