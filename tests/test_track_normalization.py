"""UCSC track normalization: strand flip, distinct-QNAME depth and empty tracks."""
import gzip
import importlib.util
import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
WORKFLOWS = ROOT / 'workflows'
SCRIPTS = ROOT / 'scripts'

HEADER = '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n'
SEQ = 'A' * 50

SAMPLE_SAM = HEADER + '\n'.join([
    '\t'.join(['r1', '0', 'chr1', '101', '60', '50M', '*', '0', '0', SEQ, '*']),
    '\t'.join(['r2', '0', 'chr1', '101', '60', '50M', '*', '0', '0', SEQ, '*']),
    '\t'.join(['r3', '16', 'chr1', '201', '60', '50M', '*', '0', '0', SEQ, '*']),
    '\t'.join(['r4', '0', 'chr1', '301', '60', '50M100N50M', '*', '0', '0', 'A' * 100, '*']),
    '\t'.join(['r5', '0', 'chr1', '501', '60', '50M', '*', '0', '0', SEQ, '*']),
    '\t'.join(['r5', '2048', 'chr1', '601', '60', '50M', '*', '0', '0', SEQ, '*']),
    '\t'.join(['r6', '99', 'chr1', '701', '60', '50M', '=', '801', '150', SEQ, '*']),
    '\t'.join(['r7', '147', 'chr1', '801', '60', '50M', '=', '701', '-150', SEQ, '*']),
    '\t'.join(['r8', '4', '*', '0', '0', '*', '*', '0', '0', SEQ, '*']),
]) + '\n'

UNIQUE_SAM = HEADER + '\n'.join([
    '\t'.join(['r1', '0', 'chr1', '101', '60', '50M', '*', '0', '0', SEQ, '*']),
    '\t'.join(['r3', '16', 'chr1', '201', '60', '50M', '*', '0', '0', SEQ, '*']),
]) + '\n'

EMPTY_SAM = HEADER

PLUS = [(100, 150, 2), (300, 350, 1), (450, 550, 1), (600, 650, 1), (700, 750, 1)]
MINUS = [(200, 250, 1), (800, 850, 1)]
FR_PAIRED_PLUS = [(100, 150, 2), (300, 350, 1), (450, 550, 1), (600, 650, 1), (700, 750, 1), (800, 850, 1)]
FR_PAIRED_MINUS = [(200, 250, 1)]
RF_PAIRED_PLUS = [(100, 150, 2), (300, 350, 1), (450, 550, 1), (600, 650, 1)]
RF_PAIRED_MINUS = [(200, 250, 1), (700, 750, 1), (800, 850, 1)]
UNIQUE_PLUS = [(100, 150, 1)]
UNIQUE_MINUS = [(200, 250, 1)]


def _require(*tools):
    missing = [t for t in tools if not shutil.which(t)]
    if missing:
        pytest.skip('missing ' + ', '.join(missing))


def _bam(path, sam):
    path.parent.mkdir(parents=True, exist_ok=True)
    samfile = path.with_suffix('.sam')
    samfile.write_text(sam)
    subprocess.run(['samtools', 'view', '-b', '-o', str(path), str(samfile)],
                   check=True, capture_output=True)
    subprocess.run(['samtools', 'index', str(path)], check=True, capture_output=True)


@pytest.fixture
def sandbox(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _bam(tmp_path / 'MAPPED/combo/sample_mapped_sorted.bam', SAMPLE_SAM)
    _bam(tmp_path / 'MAPPED/combo/sample_mapped_sorted_unique.bam', UNIQUE_SAM)
    _bam(tmp_path / 'MAPPED/combo/empty_mapped_sorted.bam', EMPTY_SAM)
    _bam(tmp_path / 'MAPPED/combo/empty_mapped_sorted_unique.bam', EMPTY_SAM)
    (tmp_path / 'genome.fa.gz').write_bytes(gzip.compress(b'>chr1\n' + b'A' * 1000 + b'\n'))
    (tmp_path / 'genome.fa').write_text('>chr1\n' + 'A' * 1000 + '\n')
    (tmp_path / 'genome.fa.fai').write_text('chr1\t1000\t6\t50\t51\n')
    (tmp_path / 'genome.chrom.sizes').write_text('chr1\t1000\n')
    (tmp_path / 'TMP').mkdir()
    for env in ('bedtools', 'perl', 'ucsc', 'samtools', 'base'):
        shutil.copy2(ROOT / 'envs' / f'{env}.yaml', tmp_path / f'{env}.yaml')
    return tmp_path


def _smk(tmp_path, stranded, paired):
    prefix = (
        "import os\n"
        "def samplecond(samples, config): return list(samples)\n"
        "def tool_params(*args): return {'OPTIONS': {}}\n"
        "def get_sortmem(w, resources): return 4\n"
        "SAMPLES = ['sample', 'empty']\n"
        "config = {}\n"
        "scombo = 'combo'\n"
        "combo = 'combo'\n"
        "rundedup = False\n"
        f"TRACKS_STRANDED = {stranded!r}\n"
        f"TRACKS_PAIRED = {paired!r}\n"
        "SETS = 'combo'\n"
        "REFERENCE = 'genome.fa.gz'\n"
        "REFDIR = '.'\n"
        f"BINS = {str(SCRIPTS)!r}\n"
        "MAXTHREAD = 1\n"
    )
    downstream = (
        "\nrule downstream:\n"
        "    input: expand(\"TRACKS/combo/{file}_mapped_{type}.{orient}.bw.trackdone\","
        " file=['sample', 'empty'], type=['sorted', 'sorted_unique'], orient=['fw', 're'])\n"
        "    output: \"downstream.done\"\n"
        "    shell: \"touch {output}\"\n"
    )
    shutil.copy2(WORKFLOWS / 'manipulate_genome.smk', tmp_path / 'manipulate_genome.smk')
    snake = tmp_path / 'Snakefile'
    snake.write_text(prefix + (WORKFLOWS / 'ucsc.smk').read_text() + downstream)
    result = subprocess.run(['snakemake', '-s', str(snake), '-j', '1', 'downstream.done'],
                            cwd=tmp_path, capture_output=True, text=True, timeout=180)
    assert result.returncode == 0, result.stdout + result.stderr
    return result


def _nf(tmp_path, stranded, paired):
    from MONSDA.Workflows import nf_paramize
    prefix = (
        'def get_always(key) {\n    return params[key]\n}\n'
        "SCOMBO = 'combo'\nCOMBO = 'combo'\nCONDITION = 'condition'\n"
        "THREADS = 1\n"
        "LONGSAMPLES = ['sample', 'empty']\n"
        f"TRACKSSTRANDED = {json.dumps(stranded)}\n"
        f"TRACKSPAIRED = {json.dumps(paired)}\n"
    )
    suffix = '''
process downstream {
    input:
    path trackdb
    output:
    path 'downstream.done'
    script:
    "touch downstream.done"
}
workflow {
    TRACKS(Channel.of('dummy'))
    downstream(TRACKS.out.trackdb)
}
'''
    workflow = tmp_path / 'test.nf'
    template = (WORKFLOWS / 'ucsc.nf').read_text().replace(
        'include { UnzipGenome; UnzipGenome_no_us } from "manipulate_genome.nf"',
        (WORKFLOWS / 'manipulate_genome.nf').read_text())
    workflow.write_text(nf_paramize(prefix + template + suffix))
    params = tmp_path / 'params.json'
    params.write_text(json.dumps(dict(
        BINS=str(SCRIPTS), TRACKSENV='ucsc', TRACKSBIN='ucsc',
        REFERENCE='genome.fa.gz', REFDIR='.', SETS='combo',
        ucsc_TRACKS_params_UCSC='', MAXTHREAD=1)))
    result = subprocess.run(['nextflow', 'run', str(workflow), '-params-file', str(params),
                             '-work-dir', str(tmp_path / 'work')],
                            cwd=tmp_path, capture_output=True, text=True, timeout=300,
                            env={**os.environ, 'JAVA_TOOL_OPTIONS': '-XX:-UseContainerSupport'})
    assert result.returncode == 0, result.stdout + result.stderr
    return result


def _bedg(path):
    with gzip.open(path, 'rt') as handle:
        return [line.rstrip('\n').split('\t') for line in handle if line.strip()]


def _check_bedg(path, expected, units):
    rows = _bedg(path)
    assert [(int(r[1]), int(r[2])) for r in rows] == [(s, e) for s, e, _ in expected]
    for row, (_, _, coverage) in zip(rows, expected):
        assert float(row[3]) == pytest.approx(coverage * 1e6 / units, rel=1e-9)


def _check_bigwig(path, expected, units):
    assert path.read_bytes()[:4] == b'\x26\xfc\x8f\x88', path
    if importlib.util.find_spec('pyBigWig') is None:
        return
    import pyBigWig
    handle = pyBigWig.open(str(path))
    try:
        intervals = handle.intervals('chr1') or ()
    finally:
        handle.close()
    assert [(s, e) for s, e, _ in intervals] == [(s, e) for s, e, _ in expected]
    for (_, _, value), (_, _, coverage) in zip(intervals, expected):
        assert value == pytest.approx(coverage * 1e6 / units, rel=1e-9)


@pytest.mark.parametrize('stranded,paired,plus,minus,uplus,uminus', [
    ('fr', 'paired', FR_PAIRED_PLUS, FR_PAIRED_MINUS, UNIQUE_PLUS, UNIQUE_MINUS),
    ('rf', 'single', MINUS, PLUS, UNIQUE_MINUS, UNIQUE_PLUS),
    ('rf', 'paired', RF_PAIRED_PLUS, RF_PAIRED_MINUS, UNIQUE_PLUS, UNIQUE_MINUS),
])
def test_snakemake_normalization(sandbox, stranded, paired, plus, minus, uplus, uminus):
    _require('snakemake', 'bedtools', 'samtools', 'bedGraphToBigWig')
    _smk(sandbox, stranded, paired)
    tracks = sandbox / 'TRACKS/combo'
    _check_bedg(tracks / 'sample_mapped_sorted.fw.norm.bedg.gz', plus, 7)
    _check_bedg(tracks / 'sample_mapped_sorted.re.norm.bedg.gz', minus, 7)
    _check_bedg(tracks / 'sample_mapped_sorted_unique.fw.norm.bedg.gz', uplus, 2)
    _check_bedg(tracks / 'sample_mapped_sorted_unique.re.norm.bedg.gz', uminus, 2)
    _check_bigwig(tracks / 'sample_mapped_sorted.fw.bw', plus, 7)
    _check_bigwig(tracks / 'sample_mapped_sorted.re.bw', minus, 7)
    _check_bigwig(tracks / 'sample_mapped_sorted_unique.fw.bw', uplus, 2)


def test_snakemake_empty_library(sandbox):
    _require('snakemake', 'bedtools', 'samtools', 'bedGraphToBigWig')
    _smk(sandbox, 'fr', 'paired')
    tracks = sandbox / 'TRACKS/combo'
    for orient in ('fw', 're'):
        assert _bedg(tracks / f'empty_mapped_sorted.{orient}.norm.bedg.gz') == []
        _check_bigwig(tracks / f'empty_mapped_sorted.{orient}.bw', [], 1)
    assert 'empty' in (sandbox / 'LOGS/combo/empty/TRACKS/ucsc/sorted_bedgtoucsc.log').read_text()


def test_snakemake_restart_with_existing_beds(sandbox):
    _require('snakemake', 'bedtools', 'samtools', 'bedGraphToBigWig')
    _smk(sandbox, 'fr', 'paired')
    bed = sandbox / 'BED/combo/sample_mapped_sorted.bed.gz'
    before = bed.stat().st_mtime_ns
    result = _smk(sandbox, 'fr', 'paired')
    assert 'Nothing to be done' in result.stdout + result.stderr
    assert bed.stat().st_mtime_ns == before


@pytest.mark.parametrize('stranded,paired,plus,minus,uplus,uminus', [
    ('fr', 'paired', FR_PAIRED_PLUS, FR_PAIRED_MINUS, UNIQUE_PLUS, UNIQUE_MINUS),
    ('rf', 'single', MINUS, PLUS, UNIQUE_MINUS, UNIQUE_PLUS),
])
def test_nextflow_normalization(sandbox, stranded, paired, plus, minus, uplus, uminus):
    _require('nextflow', 'bedtools', 'samtools', 'bedGraphToBigWig')
    _nf(sandbox, stranded, paired)
    tracks = sandbox / 'TRACKS/combo/condition'
    _check_bedg(tracks / 'sample_mapped_sorted.fw.norm.bedg.gz', plus, 7)
    _check_bedg(tracks / 'sample_mapped_sorted.re.norm.bedg.gz', minus, 7)
    _check_bedg(tracks / 'sample_mapped_sorted_unique.fw.norm.bedg.gz', uplus, 2)
    _check_bedg(tracks / 'sample_mapped_sorted_unique.re.norm.bedg.gz', uminus, 2)
    _check_bigwig(tracks / 'sample_mapped_sorted.fw.bw', plus, 7)
    _check_bigwig(tracks / 'sample_mapped_sorted.re.bw', minus, 7)
    for orient in ('fw', 're'):
        assert _bedg(tracks / f'empty_mapped_sorted.{orient}.norm.bedg.gz') == []
        _check_bigwig(tracks / f'empty_mapped_sorted.{orient}.bw', [], 1)
