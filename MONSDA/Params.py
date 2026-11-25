# Params.py ---

# Code:
# import os, sys, inspect
# # realpath() will make your script run, even if you symlink it :)
# cmd_folder = os.path.dirname(os.path.realpath(os.path.abspath( inspect.getfile( inspect.currentframe() )) ))
# if cmd_folder not in sys.path:
#     sys.path.insert(0, cmd_folder)

#
# # Use this if you want to include modules from a subfolder
# cmd_subfolder = os.path.join(os.path.dirname(os.path.realpath(os.path.abspath( inspect.getfile( inspect.currentframe() )) )),"lib")
# if cmd_subfolder not in sys.path:
#     sys.path.insert(0, cmd_subfolder)
#
# # Info:
# # cmd_folder = os.path.dirname(os.path.abspath(__file__)) # DO NOT USE __file__ !!!
# # __file__ fails if the script is called in different ways on Windows.
# # __file__ fails if someone does os.chdir() before.
# # sys.argv[0] also fails, because it doesn't not always contain the path.

import bisect
import csv
import datetime
import glob
import inspect
import itertools
import json
import os
import re
import shutil
import sys
import traceback as tb
from collections import OrderedDict, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import NamedTuple, Optional

from natsort import natsorted
from snakemake.common.configfile import load_configfile as _load_configfile

import MONSDA.Utils as mu
from MONSDA.Utils import check_run as check_run
from MONSDA.Utils import get_dict_hash, makeoutdir, setup_logger

# cmd_subfolder = [os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"../MONSDA/lib"), os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"MONSDA/lib"), os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"../lib"), os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"lib")]
# for x in cmd_subfolder:
# if x not in sys.path:
#        sys.path.insert(0, x)


try:
    scriptname = (
        os.path.basename(inspect.stack()[-1].filename)
        .replace("Run", "")
        .replace(".py", "")
    )
    log = setup_logger(scriptname)
except Exception:
    exc_type, exc_value, exc_tb = sys.exc_info()
    tbe = tb.TracebackException(
        exc_type,
        exc_value,
        exc_tb,
    )
    print("".join(tbe.format()), file=sys.stderr)


_FASTQ_SUFFIX = (
    r"(?:_L(?P<lane1>\d{1,3}))?_[Rr](?P<read>[12])"
    r"(?:_L(?P<lane2>\d{1,3}))?(?:_(?P<chunk>\d+))?\.fastq\.gz$"
)


class _FastqName(NamedTuple):
    """Parsed technical suffix of a FASTQ file name."""

    sample: str
    lane: Optional[int]
    read: Optional[int]
    chunk: Optional[int]


def _parse_fastq_name(
    filename: str, sample: Optional[str] = None
) -> Optional[_FastqName]:
    """Parse read/lane suffix from a FASTQ file name.

    Accepted examples (both orderings):
      Sample_L001_R1.fastq.gz
      Sample_R1_L001.fastq.gz
      Sample_R1_L001_001.fastq.gz

    When *sample* is given the pattern is anchored on that configured sample
    name and the plain ``Sample.fastq.gz`` stem is accepted as well.
    """
    base = os.path.basename(filename)
    if sample is not None and base == sample + ".fastq.gz":
        return _FastqName(sample, None, None, None)
    prefix = re.escape(sample) if sample is not None else r"(?P<sample>.+?)"
    match = re.match(rf"^{prefix}{_FASTQ_SUFFIX}", base)
    if not match:
        return None
    lane = match.group("lane1") or match.group("lane2")
    chunk = match.group("chunk")
    return _FastqName(
        sample if sample is not None else match.group("sample"),
        int(lane) if lane else None,
        int(match.group("read")),
        int(chunk) if chunk else None,
    )


def _is_lane_split_file(filename: str, sample: Optional[str] = None) -> bool:
    """Return True if filename matches lane split naming."""
    parsed = _parse_fastq_name(filename, sample)
    return bool(parsed and parsed.lane is not None)


def _split_index(filename: str) -> tuple:
    """Return the (lane, chunk) ordering key of a split FASTQ name."""
    parsed = _parse_fastq_name(filename)
    if not parsed:
        return (sys.maxsize, sys.maxsize)
    return (
        parsed.lane if parsed.lane is not None else sys.maxsize,
        parsed.chunk if parsed.chunk is not None else sys.maxsize,
    )


def _matches_sample_fastq(filename: str, sample: str) -> bool:
    """Match configured sample against accepted FASTQ naming schemes.

    Handles both lane orderings:
      - SAMPLE_L###_R[12] (standard Illumina)
      - SAMPLE_R[12]_L### (alternative)
    """
    return _parse_fastq_name(filename, sample) is not None


def _strip_fastq_sample_name(filename: str) -> str:
    """Strip read/lane suffix from FASTQ filename and return sample basename.

    Handles both orderings:
      - SAMPLE_L###_R[12][_###]
      - SAMPLE_R[12]_L###[_###]
    """
    parsed = _parse_fastq_name(filename)
    if parsed:
        return parsed.sample
    return re.sub(r"\.fastq\.gz$", "", os.path.basename(filename))


def _filter_lane_files_when_merged_exists(files: list) -> list:
    """Drop split or alias FASTQs if canonical merged FASTQ exists for same sample/read.

    Keeps existing behavior for non-split files and for samples without merged targets.
    Handles both lane orderings: SAMPLE_L###_R# and SAMPLE_R#_L###.
    """
    if not files:
        return files

    base_names = {os.path.basename(f) for f in files}
    keep = []

    for fp in files:
        parsed = _parse_fastq_name(fp)
        if not parsed or parsed.read is None:
            keep.append(fp)
            continue

        canonical_bn = f"{parsed.sample}_R{parsed.read}.fastq.gz"
        if os.path.basename(fp) == canonical_bn:
            keep.append(fp)
            continue
        canonical_fp = os.path.join(os.path.dirname(fp), canonical_bn)
        if canonical_bn in base_names or os.path.exists(canonical_fp):
            continue
        keep.append(fp)

    return keep


def _relative_fastq_paths(files: list) -> set:
    """Reduce FASTQ paths to paths relative to the FASTQ root."""
    return {str.join(os.sep, s.split(os.sep)[1:]) for s in files}


def _sample_ids_from_files(files: list, sample: str) -> set:
    """Reduce relative FASTQ paths to sample IDs, keeping configured names."""
    ids = set()
    for x in files:
        if _matches_sample_fastq(os.path.basename(x), os.path.basename(sample)):
            ids.add(sample)
        else:
            ids.add(
                os.path.join(
                    os.path.dirname(x),
                    _strip_fastq_sample_name(os.path.basename(x)),
                )
            )
    return ids


def _configured_samples(search: str, config: dict) -> list:
    """Return the configured sample names of one condition."""
    return natsorted(
        {
            x.replace(" ", "")
            for x in mu.get_from_dict(config["SETTINGS"], search)[0]["SAMPLES"]
        }
    )


def _fastq_dir_files(directory: str, cache: dict) -> list:
    """List the merged-preferring FASTQ files of one directory, cached per call."""
    if directory not in cache:
        try:
            names = sorted(os.listdir(directory))
        except OSError:
            names = []
        cache[directory] = _filter_lane_files_when_merged_exists(
            [os.path.join(directory, n) for n in names if n.endswith(".fastq.gz")]
        )
    return cache[directory]


def _match_configured_samples(files: list, samples: list) -> dict:
    """Map FASTQ files to the most specific configured sample they belong to."""
    bases = [os.path.basename(f) for f in files]
    matched = {}
    for _, sample in sorted(enumerate(samples), key=lambda i: (-len(i[1]), i[0])):
        for index in range(bisect.bisect_left(bases, sample), len(bases)):
            if not bases[index].startswith(sample):
                break
            if files[index] in matched:
                continue
            if _matches_sample_fastq(bases[index], sample):
                matched[files[index]] = sample
    return matched


def _discover_condition(
    search: str, config: dict, nocheck: str = None, cache: dict = None
) -> tuple:
    """Return configured samples and matched FASTQ paths of one condition."""
    if cache is None:
        cache = {}
    samples = _configured_samples(search, config)
    if nocheck is not None:
        return samples, [os.sep.join([os.sep.join(search), s]) for s in samples]
    for depth in range(len(search), len(search) - 2, -1):
        directory = os.sep.join(["FASTQ"] + list(search[0:depth]))
        files = _fastq_dir_files(directory, cache)
        matched = _match_configured_samples(files, samples)
        if matched:
            return samples, matched
        if files:
            error_msg = (
                f"No sample files found for condition {os.sep.join(search)}. "
                f"Expected to find files matching samples {samples} "
                f"in directory: {os.sep.join(['FASTQ'] + list(search))}"
            )
            log.error(scriptname + ".Params_get_samples_from_dir: " + error_msg)
            raise ValueError(error_msg)
    return samples, []


def samplesheet_to_settings(samplesheet_path: str) -> dict:
    """Read a CSV/TSV samplesheet and return a SETTINGS dict compatible with MONSDA config.

    Expected columns (case-insensitive header row):
      CONDITION  - slash-separated condition path, e.g. ``Ecoli/WT`` or ``Ecoli/WT/dummylevel``
      SAMPLE     - sample name / accession
      GROUP      - group label for differential analysis
      SEQUENCING - e.g. ``paired`` or ``single``
      REFERENCE  - path to genome FASTA (.fa, .fa.gz or .fa.bgz)
      GTF        - path to GTF annotation (.gtf, .gtf.gz or .gtf.bgz)  [optional]
      GFF        - path to GFF annotation (.gff, .gff.gz or .gff.bgz)  [optional]
      INDEX      - path to pre-built index            [optional]
      PREFIX     - mapper index prefix                [optional]
      DECOY      - path to decoy file                 [optional]
      TYPE       - sample type label                  [optional]
      BATCH      - batch label                        [optional]
      IP         - IP protocol info (for PEAKS)       [optional]

    Per-condition metadata (SEQUENCING, REFERENCE, …) only needs to be present
    on the first row for that condition; subsequent rows may leave those cells
    empty (fill-down behaviour).

    Parameters
    ----------
    samplesheet_path : str
        Absolute or relative path to the samplesheet file.

    Returns
    -------
    dict
        Nested dict suitable for assigning to ``config["SETTINGS"]``.
    """
    logid = scriptname + ".samplesheet_to_settings: "

    # --- detect delimiter ---
    with open(samplesheet_path, newline="") as fh:
        sample = fh.read(4096)

    delimiter = None
    # 1. trust the file extension
    ext = os.path.splitext(samplesheet_path)[1].lower()
    if ext in (".tsv", ".txt"):
        delimiter = "\t"
    elif ext == ".csv":
        delimiter = ","
    else:
        # 2. try the sniffer
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
            delimiter = dialect.delimiter
        except csv.Error:
            pass
        # 3. manual probe: whichever candidate appears more on the first line
        if delimiter is None:
            first_line = sample.splitlines()[0] if sample else ""
            delimiter = "\t" if first_line.count("\t") >= first_line.count(",") else ","

    settings: dict = {}

    # per-condition accumulator for fill-down metadata
    cond_meta: dict = {}

    with open(samplesheet_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        # normalise header keys to upper-case
        if reader.fieldnames is None:
            raise ValueError(
                "Samplesheet appears to be empty or has no header row: "
                + samplesheet_path
            )
        reader.fieldnames = [f.strip().upper() for f in reader.fieldnames]

        for row in reader:
            row = {
                k.strip().upper(): (v.strip() if v is not None else "")
                for k, v in row.items()
                if k
            }

            condition_str = row.get("CONDITION", "").strip()
            sample_name = row.get("SAMPLE", "").strip()
            if not condition_str or not sample_name:
                log.warning(
                    logid + "Skipping row with missing CONDITION or SAMPLE: " + str(row)
                )
                continue

            # fill-down: carry over non-empty metadata from previous rows of same condition
            if condition_str not in cond_meta:
                cond_meta[condition_str] = {}
            for key in (
                "SEQUENCING",
                "REFERENCE",
                "GTF",
                "GFF",
                "INDEX",
                "PREFIX",
                "DECOY",
                "IP",
            ):
                val = row.get(key, "")
                if val:
                    cond_meta[condition_str][key] = val

            meta = cond_meta[condition_str]

            # --- navigate / create nested dict path ---
            path_parts = [p.strip() for p in condition_str.split("/") if p.strip()]
            node = settings
            for part in path_parts:
                node = node.setdefault(part, {})

            # --- initialise leaf node on first encounter ---
            if "SAMPLES" not in node:
                node["SAMPLES"] = []
                node["GROUPS"] = []
                node["TYPES"] = []
                node["BATCHES"] = []
                node["SEQUENCING"] = meta.get("SEQUENCING", "")
                node["REFERENCE"] = meta.get("REFERENCE", "")
                node["INDEX"] = meta.get("INDEX", "")
                node["PREFIX"] = meta.get("PREFIX", "")
                node["IP"] = meta.get("IP", "")
                gtf = meta.get("GTF", "")
                gff = meta.get("GFF", "")
                node["ANNOTATION"] = {"GTF": gtf, "GFF": gff}
                decoy = meta.get("DECOY", "")
                node["DECOY"] = {decoy: ""} if decoy else {}

            node["SAMPLES"].append(sample_name)
            node["GROUPS"].append(row.get("GROUP", ""))
            node["TYPES"].append(row.get("TYPE", ""))
            node["BATCHES"].append(row.get("BATCH", ""))

    log.info(logid + "Built SETTINGS from samplesheet: " + str(list(settings.keys())))
    return settings


def inject_samplesheet_settings(configfile: str, samplesheet_path: str) -> str:
    """Load *configfile*, populate ``SETTINGS`` from *samplesheet_path* if absent,
    write the augmented config to ``<base>_with_settings.json`` and return that path.

    Parameters
    ----------
    configfile : str
        Path to the original MONSDA JSON config.
    samplesheet_path : str
        Path to the CSV/TSV samplesheet.

    Returns
    -------
    str
        Path to the written (augmented) config file.
    """
    logid = scriptname + ".inject_samplesheet_settings: "

    config = _load_configfile(configfile)

    existing = config.get("SETTINGS", {})
    # strip comment-only SETTINGS that have no SAMPLES anywhere
    has_samples = (
        any(
            isinstance(v, dict) and "SAMPLES" in v
            for cond in existing.values()
            if isinstance(cond, dict)
            for v in cond.values()
        )
        if existing
        else False
    )

    if has_samples:
        log.info(
            logid
            + "Config already contains SETTINGS with sample data; samplesheet will be ignored."
        )
        return configfile

    log.info(logid + "Populating SETTINGS from samplesheet: " + samplesheet_path)
    config["SETTINGS"] = samplesheet_to_settings(samplesheet_path)

    base, ext = os.path.splitext(configfile)
    out_path = base + "_with_settings" + (ext if ext else ".json")
    with open(out_path, "w") as fh:
        json.dump(config, fh, indent=4)
    log.info(logid + "Augmented config written to: " + out_path)
    return out_path


def _merge_lane_files(target: str, lane_candidates: list, logid: str) -> str:
    """Concatenate *lane_candidates* gzip files into *target*.

    Returns *target* on success; raises on any I/O error.
    """
    tmp = target + ".tmp"
    try:
        with open(tmp, "wb") as outfh:
            for lane_file in lane_candidates:
                with open(lane_file, "rb") as infh:
                    shutil.copyfileobj(infh, outfh)
        os.replace(tmp, target)
    except Exception:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise
    return target


def prepare_lane_split_fastqs(config: dict, max_workers: Optional[int] = None) -> int:
    """Concatenate split FASTQs into canonical files when needed.

    Paired and single-cell samples are merged into ``<id>_R1/_R2.fastq.gz``,
    single-end samples into ``<id>.fastq.gz``.

    Merges are executed in parallel using a ``ThreadPoolExecutor`` so that all
    samples/reads are processed concurrently (I/O-bound work, thread-safe).

    Parameters
    ----------
    config : dict
        Parsed MONSDA config.
    max_workers : int, optional
        Maximum number of parallel merge threads.  When ``None`` (default),
        the value is taken from ``config["MAXTHREADS"]``; if that is also
        absent the ``ThreadPoolExecutor`` built-in default is used.

    This is intentionally additive: existing canonical files are kept untouched.
    """
    logid = scriptname + ".Params_prepare_lane_split_fastqs: "
    merged_files = 0

    if max_workers is None:
        try:
            max_workers = int(config["MAXTHREADS"])
        except (KeyError, TypeError, ValueError):
            max_workers = None  # fall back to ThreadPoolExecutor default

    cache = {}
    merge_tasks: list = []
    targets = set()

    for search in [
        list(k) for k in mu.keysets_from_dict(config["SETTINGS"], "SAMPLES")
    ]:
        samples, matched = _discover_condition(search, config, None, cache)
        sequencing = str(
            mu.sub_dict(config["SETTINGS"], search).get("SEQUENCING") or ""
        ).split(",")[0]
        reads = (
            ["1", "2"]
            if "paired" in sequencing or "singlecell" in sequencing
            else ["1"]
        )
        by_sample = defaultdict(list)
        for path, sample in matched.items():
            by_sample[sample].append(path)
        for sample in samples:
            if sample not in by_sample:
                continue
            sample_dir = os.path.dirname(by_sample[sample][0])
            for read in reads:
                candidates = []
                for path in by_sample[sample]:
                    parsed = _parse_fastq_name(path, sample)
                    if parsed and parsed.read == int(read):
                        candidates.append(path)
                if not candidates:
                    continue
                if "paired" in sequencing or "singlecell" in sequencing:
                    target = os.path.join(sample_dir, f"{sample}_R{read}.fastq.gz")
                else:
                    target = os.path.join(sample_dir, f"{sample}.fastq.gz")
                if target in targets:
                    continue
                targets.add(target)
                if os.path.exists(target):
                    log.info(
                        logid
                        + f"Found split files for {sample} R{read}, but target already exists: {target} (keeping existing file)"
                    )
                    continue
                log.info(
                    logid
                    + f"Concatenating {len(candidates)} split files into {target}"
                )
                merge_tasks.append(
                    (
                        target,
                        sorted(
                            candidates,
                            key=lambda f: (_split_index(f), os.path.basename(f)),
                        ),
                    )
                )

    if not merge_tasks:
        log.debug(logid + "No lane-split FASTQ files required concatenation")
        return 0

    # --- run all merges in parallel ---
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = {
            pool.submit(_merge_lane_files, target, lanes, logid): target
            for target, lanes in merge_tasks
        }
        for future in as_completed(futures):
            target = futures[future]
            try:
                future.result()
                merged_files += 1
            except Exception:
                exc_type, exc_value, exc_tb = sys.exc_info()
                tbe = tb.TracebackException(exc_type, exc_value, exc_tb)
                log.error(
                    logid
                    + f"Failed to merge lane files into {target}: "
                    + "".join(tbe.format())
                )

    log.info(logid + f"Created {merged_files} concatenated lane-merged FASTQ files")
    return merged_files


@check_run
def get_samples(config: dict) -> list():
    """Check and return samples according to sample list on config.json

    Parameters
    ----------
    config : dict
        config json file parsed by snakemake.common.configfile

    Returns
    -------
    list
        list of samples found and checked
    """
    logid = scriptname + ".Params_get_samples: "
    cache = {}
    RETSAMPLES = list()
    for search in [
        list(k) for k in mu.keysets_from_dict(config["SETTINGS"], "SAMPLES")
    ]:
        samples, matched = _discover_condition(search, config, None, cache)
        log.debug(logid + "SAMPLECHECK: " + str(matched))
        if not matched:
            continue
        sequencing = str(
            mu.sub_dict(config["SETTINGS"], search).get("SEQUENCING") or ""
        ).split(",")[0]
        log.debug(logid + "PAIRED: " + sequencing)
        by_sample = defaultdict(list)
        for path, sample in matched.items():
            by_sample[sample].append(path)
        for sample in samples:
            if sample not in by_sample:
                continue
            f = sorted(_relative_fastq_paths(by_sample[sample]))
            if "paired" not in sequencing:
                RETSAMPLES.extend(
                    sorted(
                        _sample_ids_from_files(
                            f, os.sep.join(list(search) + [sample])
                        )
                    )
                )
                log.debug(logid + "PAIREDSAMPLES: " + str(f))
            else:
                RETSAMPLES.extend([x.replace(".fastq.gz", "") for x in f])
    log.debug(logid + "SAMPLETEST: " + str(RETSAMPLES))
    if len(RETSAMPLES) < 1:
        log.error(logid + "No samples found, please check config file")
        sys.exit(1)

    log.debug(logid + "SAMPLES: " + str(RETSAMPLES))
    return RETSAMPLES


@check_run
def get_samples_postprocess(config: dict, subwork: str) -> list:
    """Check and return samples according to sample list on config.json

    Parameters
    ----------
    config : dict
        config json file parsed by snakemake.common.configfile
    subwork : str
        subworkflow to check

    Returns
    -------
    list
        list of samples found and checked
    """
    logid = scriptname + ".Params_get_samples_postprocess: "
    cache = {}
    RETSAMPLES = list()
    for search in [
        list(k) for k in mu.keysets_from_dict(config["SETTINGS"], "SAMPLES")
    ]:
        samples, matched = _discover_condition(search, config, None, cache)
        by_sample = defaultdict(list)
        for path, sample in matched.items():
            by_sample[sample].append(path)
        for sample in samples:
            qualified = os.sep.join(list(search) + [sample])
            if not len(
                mu.get_from_dict(config[subwork], conditiononly(qualified, config))
            ):
                continue
            if config[subwork].get("EXCLUDE") and sample in config[subwork]["EXCLUDE"]:
                continue
            if sample not in by_sample:
                continue
            f = sorted(_relative_fastq_paths(by_sample[sample]))
            log.debug(logid + "SAMPLECHECK: " + str(f))
            RETSAMPLES.extend(sorted(_sample_ids_from_files(f, qualified)))
    log.debug(logid + "SAMPLETEST: " + str(RETSAMPLES))
    if len(RETSAMPLES) < 1:
        log.error(
            logid
            + "No samples found for "
            + str(subwork)
            + ", please check config file"
        )
        sys.exit(1)

    log.debug(logid + "SAMPLES: " + str(RETSAMPLES))
    return RETSAMPLES


@check_run
def check_samples(config: dict) -> bool:
    """Check if all samples are available in FASTQ folder

    Parameters
    ----------
    config : dict
        config json file parsed by snakemake.common.configfile

    Returns
    -------
    bool
    """
    logid = scriptname + ".Params_check_samples: "
    cache = {}
    for search in [
        list(k) for k in mu.keysets_from_dict(config["SETTINGS"], "SAMPLES")
    ]:
        samples = _configured_samples(search, config)
        directory = os.sep.join(["FASTQ"] + list(search))
        log.debug(logid + "SEARCHING: " + directory)
        matched = _match_configured_samples(_fastq_dir_files(directory, cache), samples)
        log.debug(logid + "SAMPLECHECK: " + str(matched))
        found = set(matched.values())
        if any(sample not in found for sample in samples):
            return False
    return True


@check_run
def download_samples(config: dict) -> list:
    """Download samples from config file

    Parameters
    ----------
    config : dict
        config json file parsed by snakemake.common.configfile

    Returns
    -------
    list
        list of samples found and checked
    """

    logid = scriptname + ".Params_download_samples: "
    SAMPLES = [os.path.join(x) for x in sampleslong(config, nocheck=True)]
    log.debug(logid + "DOWNLOAD_SAMPLES_LONG: " + str(SAMPLES))
    return SAMPLES


BASECALL_SIGNAL_SUFFIX = {"guppy": ".fast5", "dorado": ".pod5"}


def basecall_caller(config: dict) -> Optional[str]:
    """Return the basecaller env the config pins, None if it does not pin one."""
    basecall = config.get("BASECALL")
    basecall = basecall if isinstance(basecall, dict) else {}
    tools = [x for x in (basecall.get("TOOLS") or {}) if x in BASECALL_SIGNAL_SUFFIX]
    for env in (config.get("BASECALLENV"), basecall.get("ENV")):
        if env in BASECALL_SIGNAL_SUFFIX:
            return env
    return tools[0] if len(tools) == 1 else None


def _basecall_checked(path: str, sample: str, suffix: str) -> str:
    """Validate one resolved raw signal file and return its path."""
    where = "raw signal input " + path + " for sample " + sample + " "
    if not os.path.isfile(path):
        raise ValueError(where + "is not a regular file")
    if not path.endswith(suffix):
        raise ValueError(where + "does not end in " + suffix)
    if os.path.getsize(path) == 0:
        raise ValueError(where + "is empty")
    return path


def _basecall_list(source: str, sample: str, suffix: str, rawdir: str) -> list:
    """Resolve the loose raw signal files listed in an explicit sample list."""
    with open(source, "r") as fh:
        entries = [x.strip() for x in fh if x.strip()]
    files = list()
    for number, name in enumerate(entries, 1):
        where = "entry " + str(number) + " of " + source + " for sample " + sample
        if os.path.isabs(name):
            raise ValueError(where + " is absolute, list paths relative to " + rawdir)
        path = os.path.normpath(os.path.join(rawdir, name))
        if os.path.relpath(path, rawdir).startswith(os.pardir):
            raise ValueError(where + " points outside " + rawdir + ": " + name)
        if path in files:
            raise ValueError("duplicate entry " + name + " in " + source + " for sample " + sample)
        files.append(_basecall_checked(path, sample, suffix))
    if not files:
        raise ValueError("list " + source + " for sample " + sample + " names no raw signal file")
    return files


def _basecall_files(sample: str, suffix: str, strict: bool = True) -> Optional[list]:
    """Raw signal files of one configured sample, resolved from its single source."""
    condition, name = os.path.split(sample)
    rawdir = os.path.join("RAW", condition)
    sources = [
        ("file", os.path.join(rawdir, name + suffix)),
        ("directory", os.path.join(rawdir, name)),
        ("list", os.path.join(rawdir, name + suffix + ".list")),
    ]
    found = [
        (kind, path)
        for kind, path in sources
        if (os.path.isdir(path) if kind == "directory" else os.path.isfile(path))
    ]
    if not found:
        if not strict:
            return None
        raise ValueError(
            "no raw signal source for sample " + sample + ", expected exactly one of "
            + ", ".join(path for _, path in sources)
        )
    if len(found) > 1 and strict:
        raise ValueError(
            "ambiguous raw signal source for sample " + sample + ": "
            + ", ".join(path for _, path in found)
            + " all exist, keep exactly one of "
            + ", ".join(path for _, path in sources)
        )
    kind, source = found[0]
    if kind == "file":
        return [_basecall_checked(source, sample, suffix)]
    if kind == "list":
        return _basecall_list(source, sample, suffix, rawdir)
    files = [
        _basecall_checked(os.path.join(root, entry), sample, suffix)
        for root, _dirs, names in os.walk(source)
        for entry in names
        if entry.endswith(suffix)
    ]
    if not files:
        raise ValueError(
            "raw signal directory " + source + " for sample " + sample
            + " contains no " + suffix + " files"
        )
    return files


def basecall_signal_map(config: dict, caller: str = None) -> dict:
    """Map condition qualified sample IDs to their raw signal files."""
    caller = caller or basecall_caller(config)
    if caller not in BASECALL_SIGNAL_SUFFIX:
        raise ValueError(
            "no basecaller selected, set BASECALLENV or a single BASECALL TOOLS entry to one of "
            + ", ".join(sorted(BASECALL_SIGNAL_SUFFIX))
        )
    suffix = BASECALL_SIGNAL_SUFFIX[caller]
    samples = natsorted(os.path.join(x) for x in sampleslong(config, nocheck=True))
    if not samples:
        raise ValueError("no samples configured for BASECALL")
    signalmap = OrderedDict((x, natsorted(_basecall_files(x, suffix))) for x in samples)
    owner = dict()
    for sample, files in signalmap.items():
        for path in files:
            real = os.path.realpath(path)
            if owner.setdefault(real, sample) != sample:
                raise ValueError(
                    "raw signal input " + path + " is claimed by samples "
                    + owner[real] + " and " + sample
                )
    return dict(signalmap)


@check_run
def basecall_samples(config: dict) -> list:
    """Check and return samples according to sample list on config.json

    Parameters
    ----------
    config : dict
        config json file parsed by snakemake.common.configfile

    Returns
    -------
    list
        list of samples found and checked
    """
    logid = scriptname + ".Params_basecall_samples: "
    caller = basecall_caller(config)
    try:
        RETSAMPLES = (
            list(basecall_signal_map(config, caller))
            if caller
            else [
                sample
                for sample in (
                    os.path.join(x) for x in sampleslong(config, nocheck=True)
                )
                if any(
                    _basecall_files(sample, suffix, strict=False)
                    for suffix in sorted(set(BASECALL_SIGNAL_SUFFIX.values()))
                )
            ]
        )
    except ValueError as err:
        log.error(logid + str(err))
        sys.exit(1)
    log.debug(logid + "SAMPLETEST: " + str(RETSAMPLES))
    if len(RETSAMPLES) < 1:
        log.error(logid + "No samples found, please check config file")
        sys.exit(1)

    log.debug(logid + "SAMPLES: " + str(RETSAMPLES))
    return natsorted(RETSAMPLES)


@check_run
def get_conditions(config: dict, stage: str = "SETTINGS") -> list:
    """Get conditions from config file based on stage

    Parameters
    ----------
    config : dict
        config json file parsed by snakemake.common.configfile
    stage : str
        stage to check

    Returns
    -------
    list
        list of conditions configured in config file
    """

    logid = scriptname + ".Params_conditions: "
    ret = list()
    for k in mu.keysets_from_dict(config[stage], "SAMPLES"):
        ret.append(k)
    log.debug(logid + str(ret))
    return sorted(list(set(ret)))


@check_run
def get_samples_from_dir(search: str, config: dict, nocheck: str = None) -> list:  # type: ignore
    """Get samples from directory based on search string

    Parameters
    ----------
    search : str
        search string
    config : dict
        config json file parsed by snakemake.common.configfile
    nocheck : str, default=None
        nocheck string

    Returns
    -------
    list
        list of samples found
    """
    logid = scriptname + ".Params_get_samples_from_dir: "
    samples, matched = _discover_condition(search, config, nocheck, {})
    log.debug(logid + f"Samples: {samples}, Search: {search}, Check: {nocheck}")
    if nocheck is not None:
        return list(set(matched))
    if matched:
        sequencing = str(
            mu.sub_dict(config["SETTINGS"], search).get("SEQUENCING") or ""
        ).split(",")[0]
        paths = natsorted(matched, key=lambda y: y.lower())
        if "paired" in sequencing or "singlecell" in sequencing:
            ret = list(
                dict.fromkeys(
                    os.sep.join(
                        [
                            os.sep.join(os.path.dirname(s).split(os.sep)[1:]),
                            _strip_fastq_sample_name(os.path.basename(s)),
                        ]
                    )
                    for s in paths
                )
            )
        else:
            ret = list(
                dict.fromkeys(
                    os.sep.join(s.split(os.sep)[1:]).replace(".fastq.gz", "")
                    for s in paths
                )
            )
        log.debug(logid + "RETURN: " + str(ret))
        return ret
    log.error(logid + "NO SAMPLES FOUND")
    return list()


@check_run
def sampleslong(config: dict, nocheck: str = None) -> list:  # type: ignore
    """Get samples from config file with condition and sample name

    Parameters
    ----------
    config : dict, optional
        config json file parsed by snakemake.common.configfile

    Returns
    -------
    list
        list of samples found
    """
    logid = scriptname + ".Params_sampleslong: "
    log.debug(logid + f"Check: {nocheck}")
    tosearch = list()
    ret = list()
    for k in mu.keysets_from_dict(config["SETTINGS"], "SAMPLES"):
        tosearch.append(list(k))
    log.debug(logid + "tosearch: " + str(tosearch))
    for search in tosearch:
        ret.extend(get_samples_from_dir(search, config, nocheck))
    log.debug(logid + str(ret))
    return ret


@check_run
def get_placeholder(config: dict) -> list:
    """Get placeholder from config file

    Parameters
    ----------
    config : dict
        config json file parsed by snakemake.common.configfile

    Returns
    -------
    list
        list of placeholders found
    """
    ret = list()
    if "PH" in (config):
        for x in config["PH"]:
            ret.append(str(x))
    else:
        ret.append("_")
    return ret


@check_run
def get_cutoff_as_string(config: dict, subwork: str, cf: str) -> str:
    """Get cutoff as string from config file

    Parameters
    ----------
    config : dict
        config json file parsed by snakemake.common.configfile
    subwork : str
        workflow part
    cf : str
        cutoff, either 'padj' or 'lfc'

    Returns
    -------
    str
        configured cutoffs or defaults
    """
    logid = scriptname + ".get_cutoff: "
    defaults = {"padj": ".05", "lfc": "1.5"}
    deprecated = {"padj": ["pvalue", "pval"]}
    cutoffs = config[subwork].get("CUTOFFS") or dict()
    cutoff = cutoffs.get(cf)
    if cutoff is None:
        for alias in deprecated.get(cf, []):
            if cutoffs.get(alias) is not None:
                cutoff = cutoffs[alias]
                log.warning(
                    logid
                    + "CUTOFFS key '"
                    + alias
                    + "' in "
                    + str(subwork)
                    + " is deprecated, cutoffs are applied to adjusted p-values, please rename it to '"
                    + cf
                    + "'"
                )
                break
    if cutoff is None:
        if cf not in defaults:
            log.error(logid + "No default cutoff known for " + str(cf))
        cutoff = defaults.get(cf, "1.5")
    cutoff = str(cutoff)
    log.info(logid + "CUTOFFS: " + str(cf) + ":" + cutoff)
    return cutoff


@check_run
def get_summary_dirs(config: dict) -> list:
    """Get summary directories from config file

    Parameters
    ----------
    config : dict
        config json file parsed by snakemake.common.configfile

    Returns
    -------
    list
        list of summary directories
    """
    logid = scriptname + ".get_summary_dirs: "
    ret = list()
    for work, tools in config["WORKFLOWS"].items():
        for tool in tools:
            ret.append(f"{work}/{tool.upper()}")
    log.debug(logid + str(ret))
    return ret


@check_run
def get_summary_files(config: dict) -> list:
    """Get summary files from config file

    Parameters
    ----------
    config : dict
        config json file parsed by snakemake.common.configfile

    Returns
    -------
    list
        list of summary files
    """
    logid = scriptname + ".get_summary_files: "
    ret = list()
    for work, tools in config["WORKFLOWS"].items():
        for tool in tools:
            log.info(logid + "make summary of " + str(work) + " - " + str(tool))
            for f in glob.glob(f"{work}/{tool.upper()}/Sig*"):
                # for f in glob.glob(f"{work}/{tool.upper()}/*"):
                ret.append(f)
    log.debug(logid + str(ret))
    return ret


@check_run
def create_skeleton(runner: str, skeleton: str = None) -> None:
    """Create skeleton directories

    Parameters
    ----------
    runner : str
        workflow part
    skeleton : str, optional
        set to create skeleton directories, by default None
    """
    logid = scriptname + ".Params_create_skeleton: "
    if skeleton:
        for subdir in ["SubSnakes", "RAW", "FASTQ", "LOGS", "TMP", "JOB"]:
            mu.makeoutdir(subdir)
            sys.exit(
                "Skeleton directories created, please add files and rerun without --skeleton option"
            )
    else:
        for subdir in [runner, "LOGS", "TMP", "JOBS"]:
            mu.makeoutdir(subdir)
        if os.path.isfile(os.path.abspath("JOBS" + os.sep + scriptname + ".commands")):
            ts = str(
                datetime.datetime.fromtimestamp(
                    os.path.getmtime(
                        os.path.abspath("JOBS" + os.sep + scriptname + ".commands")
                    )
                ).strftime("%Y%m%d_%H_%M_%S")
            )
            shutil.copy2(
                "JOBS" + os.sep + scriptname + ".commands",
                "JOBS" + os.sep + scriptname + "_" + ts + ".commands",
            )
            open("JOBS" + os.sep + scriptname + ".commands", "w").close()


COUNT_STAGES = ("COUNTING", "DE", "DEU", "DAS", "DTU")


def normalize_count_options(options: dict) -> dict:
    """Normalize the public OPTIONS.QUANT key, folding the legacy OPTIONS.COUNT into it

    Parameters
    ----------
    options : dict
        tool OPTIONS dict, not modified

    Returns
    -------
    dict
        copy of options with the public QUANT key
    """
    logid = scriptname + ".Params_normalize_count_options: "
    if not mu.dict_inst(options):
        return options
    opts = OrderedDict(options)
    quant_set = bool(str(opts.get("QUANT") or "").strip())
    count_set = bool(str(opts.get("COUNT") or "").strip())
    if count_set and not quant_set:
        opts["QUANT"] = opts["COUNT"]
    elif not quant_set and "QUANT" not in opts and "COUNT" in opts:
        opts["QUANT"] = opts["COUNT"]
    opts.pop("COUNT", None)
    if "COUNT" in options and count_set != quant_set:
        used = "OPTIONS.QUANT" if quant_set else "OPTIONS.COUNT"
        log.warning(
            logid
            + "OPTIONS.QUANT is the documented key, OPTIONS.COUNT is the legacy fallback; using "
            + used
        )
    return opts


@check_run
def tool_params(sample: str, runstate: str, config: dict, subconf: dict, tool: str = None) -> dict:  # type: ignore
    """Get tool parameters from config file for a specific sample

    Parameters
    ----------
    sample : str
        sample name
    runstate : str
        part of workflow
    config : dict
        config json file parsed by snakemake.common.configfile
    subconf : dict
        subconfig
    tool : str, optional
        tool name, by default None

    Returns
    -------
    dict
        tool parameters
    """
    logid = scriptname + ".Params_tool_params: "
    log.debug(logid + "Samples: " + str(sample))
    mp = OrderedDict()
    x = sample.split(os.sep)[:-1]
    if runstate is None:
        runstate = runstate_from_sample([sample], config)[0]
    if runstate not in x:
        x.append(runstate)
    log.debug(logid + str([sample, runstate, subconf, x]))
    if "_" in tool:
        tool = tool.split("_")[0]
    mp = (
        mu.sub_dict(config[subconf], x)[tool]
        if tool
        else mu.sub_dict(config[subconf], x)
    )
    if tool and subconf in COUNT_STAGES and mu.dict_inst(mp) and "OPTIONS" in mp:
        mp = OrderedDict(mp)
        mp["OPTIONS"] = normalize_count_options(mp["OPTIONS"])
    log.debug(logid + "DONE: " + str(mp))
    return mp


@check_run
def setting_per_sample(sample: str, runstate: str, config: dict, setting: str, subconf: str = None) -> str:  # type: ignore
    """Set parameters per sample

    Parameters
    ----------
    sample : str
        sample name
    runstate : str
        condition tree part
    config : dict
        config json file parsed by snakemake.common.configfile
    setting : str
        workflow part
    subconf : str, optional
        subconfig, by default None

    Returns
    -------
    str
        setting per sample
    """
    logid = scriptname + ".Params_setting_per_sample: "
    log.debug(logid + "Samples: " + str(sample))
    sets = None
    x = sample.split(os.sep)[2:-1]
    if runstate is None:
        runstate = runstate_from_sample([sample], config)[0]
    if runstate not in x:
        x.append(runstate)
    subsetting = mu.sub_dict(config["SETTINGS"], x).get(setting)

    if setting == "ANNOTATION":  # Special case is annotation
        subsetting = subsetting.get(
            "GTF", subsetting.get("GFF")
        )  # by default GTF format will be used

    if subconf:  # check specific setting for workflow part
        subset = (
            config[subconf].get(setting)
            if config[subconf].get(setting)
            else mu.sub_dict(config[subconf], x).get(setting)
        )

    # Define which final setting is returned
    sets = subset if subset else setting

    return sets


@check_run
def get_reps(samples: list, config: dict, analysis: str, process: str = "smk") -> str:
    """get reps from samples

    Parameters
    ----------
    samples : list
        samples
    config : dict
        config json file parsed by snakemake.common.configfile
    analysis : str
        workflow part
    process : str, optional
        snakemake or nextflow process, by default "smk"

    Returns
    -------
    str
        condition string for differential analysis per sample
    """
    logid = scriptname + ".Params_get_reps: "
    log.debug(logid + "Samples: " + str(samples))
    ret = defaultdict(list)
    for sample in samples:
        if process == "smk":
            checkid = sample.split(os.sep)[4]
            scond = (
                sample.split(os.sep)[4:-1]
                if config["SETTINGS"].get(checkid, False)
                else sample.split(os.sep)[3:-1]
            )
        else:
            scond = sample.split(os.sep)[:-1]
        log.debug(logid + "WORKING ON: " + str(sample) + " CONDITION: " + str(scond))
        partconf = mu.sub_dict(config["SETTINGS"], scond)
        log.debug(logid + "CONF: " + str(partconf))

        Ex = config[analysis].get("EXCLUDE")
        if Ex and sample.split(os.sep)[-1] in Ex:
            continue
        if process == "smk":
            ret["reps"].append(sample)
        else:
            ret["reps"].append(os.path.basename(sample))
        wcfile = (
            sample.split(os.sep)[-1]
            .replace("_mapped_sorted_unique.counts.gz", "")
            .replace("_mapped_sorted_unique_dedup.counts.gz", "")
        )
        idx = partconf["SAMPLES"].index(wcfile)
        scond.append(sample.split(os.sep)[-1])
        ret["pairs"].append(checkpaired_rep([str.join(os.sep, scond)], config))
        ret["conds"].append(partconf["GROUPS"][idx])
        if (
            "BATCHES" in partconf
            and len(partconf["BATCHES"]) >= idx + 1
            and str(partconf["BATCHES"][idx]) != ""
        ):
            ret["batches"].append(str(partconf["BATCHES"][idx]).replace(",", "_"))
        else:
            ret["batches"].append("1")
        if (
            "TYPES" in partconf
            and len(partconf["TYPES"]) >= idx + 1
            and str(partconf["TYPES"][idx]) != ""
        ):
            ret["types"].append(str(partconf["TYPES"][idx]).replace(",", "_"))
        else:
            ret["types"].append("std")

    rets = "-r " + str.join(",", ret["reps"])
    rets += " -c " + str.join(",", ret["conds"])
    rets += " -t " + str.join(",", ret["types"])
    rets += " -b " + str.join(",", ret["batches"])
    rets += " --paired " + str.join(",", ret["pairs"]) if "pairs" in ret else ""

    log.debug(logid + "RETURN: " + str(rets))
    return rets


@check_run
def get_reps_file(
    samples: list, config: dict, analysis: str, process: str = "smk"
) -> str:
    """get reps from samples and write to file

    Parameters
    ----------
    samples : list
        samples
    config : dict
        config
    analysis : str
        workflow part
    process : str, optional
        , by default "smk"
        snakemake or nextflow process
    Returns
    -------
    str
        condition string for differential analysis per sample in file
    """
    argstr = get_reps(samples, config, analysis, process)
    toks = argstr.split()
    lines = []
    for i in range(0, len(toks), 2):
        lines.append(toks[i] + "=" + toks[i + 1])
    content = "\n".join(lines) + "\n"
    fn = os.path.join(
        "TMP", "reps_" + analysis + "_" + get_dict_hash({"c": content}) + ".args"
    )
    if not os.path.exists(fn):
        makeoutdir("TMP")
        with open(fn + "." + str(os.getpid()), "w") as fh:
            fh.write(content)
        os.replace(fn + "." + str(os.getpid()), fn)
    return "@" + fn


@check_run
def get_diego_samples(samples: list, config: dict, analysis: str) -> str:
    """samples for diego analysis

    Parameters
    ----------
    samples : list
        samples
    config : dict
        config json file parsed by snakemake.common.configfile
    analysis : str
        workflow part

    Returns
    -------
    str
        samples for diego analysis
    """
    logid = scriptname + ".Params_get_diego_samples: "
    log.debug(logid + "Samples: " + str(samples))
    ret = defaultdict(list)
    for sample in samples:
        scond = sample.split(os.sep)[4:-1]
        log.debug(logid + "WORKING ON: " + str(sample) + " CONDITION: " + str(scond))
        wcfile = str.join("-", sample.split(os.sep)[-4:]).replace(
            "_mapped_sorted_unique.counts.gz", ""
        )
        ret[wcfile].append(sample)

    log.debug(logid + "RETURN: " + str(ret))

    slist = ""
    for key, val in ret.items():
        slist += str(key) + "\t"
        slist += "\t".join(val)
        slist += os.linesep

        log.debug(logid + "RETURN: " + str(slist))
    return slist


@check_run
def get_diego_groups(samples: list, config: dict, analysis: str) -> str:
    """get diego groups

    Parameters
    ----------
    samples : list
        samples
    config : dict
        config json file parsed by snakemake.common.configfile
    analysis : str
        workflow part

    Returns
    -------
    str
        groups for diego analysis
    """
    logid = scriptname + ".Params_get_diego_groups: "
    log.debug(logid + "Samples: " + str(samples))
    ret = defaultdict(list)
    for sample in samples:
        scond = sample.split(os.sep)[4:-1]
        log.debug(logid + "WORKING ON: " + str(sample) + " CONDITION: " + str(scond))
        partconf = mu.sub_dict(config["SETTINGS"], scond)
        log.debug(logid + "CONF: " + str(partconf))
        wcfile = str.join("-", sample.split(os.sep)[-4:]).replace(
            "_mapped_sorted_unique.counts.gz", ""
        )
        checkfile = sample.split(os.sep)[-1].replace(
            "_mapped_sorted_unique.counts.gz", ""
        )
        idx = partconf["SAMPLES"].index(checkfile)
        cond = partconf["GROUPS"][idx]
        ret[cond].append(wcfile)

    slist = ""
    for key, val in ret.items():
        slist += str(key) + "\t"
        slist += "\t".join(val)
        slist += os.linesep
    log.debug(logid + "RETURN: " + str(slist))
    return slist


@check_run
def env_bin_from_config2(samples: list, config: dict, subconf: str) -> tuple:
    """get env and bin from config

    Parameters
    ----------
    samples : list
        samples
    config : dict
        config json file parsed by snakemake.common.configfile
    subconf : str
        part of config

    Returns
    -------
    tuple
        env and bin for workflow part and tool
    """
    logid = scriptname + ".Params_env_bin_from_config2: "
    for s in samples:
        log.debug(logid + "S: " + str(s))
        log.debug(logid + "C: " + str(conditiononly(s, config)))
        check = conditiononly(s, config)

        for k in mu.get_from_dict(config[subconf], check):
            if "BIN" in k:
                mb = k["BIN"]
            else:
                mb = ""
            if "ENV" in k:
                me = k["ENV"]
            else:
                me = ""
        log.debug(logid + str([str(mb), str(me)]))
    return mb, me


@check_run
def env_bin_from_config(config: dict, subconf: str) -> tuple:
    """env and bin from config

    Parameters
    ----------
    config : dict
        config json file parsed by snakemake.common.configfile
    subconf : str
        workflow part

    Returns
    -------
    tuple
        env and bin for workflow part and tool
    """
    logid = scriptname + ".Params_env_bin_from_config: "
    envkey = subconf + "ENV"
    binkey = subconf + "BIN"
    me = config[envkey]
    mb = config[binkey]
    log.debug(logid + str([str(mb), str(me)]))
    return mb, me


@check_run
def sample_from_path(path: str) -> str:
    """sample from path

    Parameters
    ----------
    path : str
        path to sample in workdir

    Returns
    -------
    str
        sample name
    """
    ret = str(os.path.join(os.path.split(str(path))[-1]))
    return ret


def sample_from_wildcard(sample: str, samples: list) -> str:
    """Resolve a wildcard sample identifier to the matching entry of samples.

    Parameters
    ----------
    sample : str
        sample or condition qualified sample name taken from a wildcard
    samples : list
        sample names as returned by sampleslong

    Returns
    -------
    str
        the entry of samples that names exactly this sample
    """
    wanted = str(sample)
    candidates = list(dict.fromkeys(str(x) for x in samples))
    if os.sep in wanted:
        if wanted in candidates:
            return wanted
        raise ValueError(
            "sample " + wanted + " not found in configured samples " + str(samples)
        )
    matches = [x for x in candidates if os.path.basename(x) == wanted]
    if not matches:
        raise ValueError(
            "sample " + wanted + " not found in configured samples " + str(samples)
        )
    if len(matches) > 1:
        raise ValueError(
            "sample "
            + wanted
            + " is ambiguous, candidates: "
            + ", ".join(matches)
            + "; use a condition qualified sample ID"
        )
    return matches[0]


@check_run
def runstate_from_sample(sample: list, config: dict) -> list:
    """all runstates for sample in config

    Parameters
    ----------
    sample : str
        sample name
    config : dict
        config json file parsed by snakemake.common.configfile

    Returns
    -------
    list
        all runstates for sample
    """
    logid = scriptname + ".Params_runstate_from_sample: "
    ret = list()
    for s in sample:
        if len(s.split(os.sep)) < 2:
            s = samplecond(s, config)
        if len(mu.get_from_dict(config["SETTINGS"], s.split(os.sep))) < 1:
            s = os.path.dirname(s)
            n = s.split(os.sep)[-1]
        log.debug(logid + "SAMPLE: " + s)
        try:
            c = mu.get_from_dict(config["SETTINGS"], s.split(os.sep))[0]
        except (KeyError, IndexError, TypeError):
            c = None
        log.debug(logid + "SETTINGS: " + str(c))
        if mu.dict_inst(c):
            if not c.get("SAMPLES"):
                for k, v in c.items():
                    log.debug(
                        logid + "k: " + str(k) + ", v: " + str(v) + " c: " + str(c)
                    )
                    if mu.dict_inst(v) and v.get("SAMPLES"):
                        if k not in ret:
                            ret.append(k)
            else:
                log.debug(logid + "c: " + str(c))
                ret.extend(s.split(os.sep))
        else:
            log.debug(logid + "c: " + str(c) + ", k: " + str(s.split(os.sep)[-1]))
            k = s.split(os.sep)[-1]
            ret.append(k)
    log.debug(logid + "RETURN: " + str(ret))
    return ret


@check_run
def samplecond(sample: str, config: dict) -> list:
    """condition for sample based on config

    Parameters
    ----------
    sample : str
        sample name
    config : dict
        config json file parsed by snakemake.common.configfile

    Returns
    -------
    list
        condition for sample
    """
    logid = scriptname + ".Params_samplecond: "
    ret = list()
    for s in sample:
        log.debug(logid + str(s))
        s = s.replace(".fastq.gz", "")
        check = s.split(os.sep)
        checkdir = check[:-1]
        sname = check[-1]
        tmplist = checkdir
        log.debug(logid + "CHECK: " + str(checkdir))
        for r in runstate_from_sample([s], config):
            if r not in tmplist:
                tmp = check[:-1]
                tmp.append(r)
                if sname in mu.get_from_dict(config["SETTINGS"], tmp)[0].get("SAMPLES"):
                    tmplist.append(r)
        log.debug(logid + "TMPLIST: " + str(tmplist))
        paired = checkpaired([s], config)
        if any(p in paired for p in ["paired", "singlecell"]):
            s = re.sub(r"_[r|R|][1|2]\.", "", s)
        r = os.sep.join(tmplist)
        if r not in s:
            ret.append(os.sep.join([r, os.path.basename(s)]))
        else:
            ret.append(os.sep.join([os.path.dirname(s), os.path.basename(s)]))
    log.debug(logid + "RETURN: " + str(ret))
    return ret


@check_run
def conditiononly(sample: str, config: dict) -> list:
    """only condition for sample

    Parameters
    ----------
    sample : str
        sample name
    config : dict
        config json file parsed by snakemake.common.configfile

    Returns
    -------
    list
        only condition for sample
    """
    logid = scriptname + ".Params_conditiononly: "
    ret = list()
    check = sample.split(os.sep)
    checkdir = check[:-1]
    sname = check[-1]
    ret.extend(checkdir)
    log.debug(logid + "CHECK: " + str(checkdir))
    for r in runstate_from_sample([sample], config):
        log.debug(logid + "runstate " + str(r))
        if r not in ret:
            tmp = list()
            tmp.extend(
                ret
            )  # this will take only the first occurence of sample in settings, should anyways never happen to have the same sample in different subsettings with differing pairedness
            tmp.append(r)
            log.debug(logid + "tmp: " + str(tmp))
            if len(
                mu.get_from_dict(config["SETTINGS"], tmp)
            ) > 0 and sname in mu.get_from_dict(config["SETTINGS"], tmp)[0].get(
                "SAMPLES"
            ):
                ret.append(r)
    log.debug(logid + "ret: " + str(ret))
    return ret


@check_run
def checkpaired(sample: list, config: dict) -> str:
    """check if sample is paired

    Parameters
    ----------
    sample : list
        samples
    config : dict
        config json file parsed by snakemake.common.configfile

    Returns
    -------
    str
        pairedness of sample
    """
    logid = scriptname + ".Params_checkpaired: "
    paired = ""
    for s in sample:
        log.debug(logid + "SAMPLE: " + str(s))
        check = conditiononly(s, config)
        log.debug(logid + "CHECK: " + str(check))
        p = mu.sub_dict(config["SETTINGS"], check)
        if p:
            paired = p.get("SEQUENCING")
            paired = paired.split(",")[0] if "," in paired else paired
        else:
            return None
        # Per sample paired, not implemented yet
        # pairedlist = p.get('SEQUENCING')
        # samplelist = p.get('SAMPLES')
        # x = samplelist.index(s.split(os.sep)[-1])
        # paired = pairedlist[x]
    log.debug(logid + "SEQUENCING: " + str(paired))
    return paired


@check_run
def checkpaired_rep(sample: list, config: dict) -> str:
    """check replicates for pairedness

    Parameters
    ----------
    sample : list
        samples
    config : dict
        config json file parsed by snakemake.common.configfile

    Returns
    -------
    str
        pairedness of replicates
    """
    logid = scriptname + ".Params_checkpaired_rep: "
    log.debug(logid + "SAMPLE: " + str(sample))
    ret = list()
    for s in sample:
        check = conditiononly(s, config)
        p = mu.sub_dict(config["SETTINGS"], check)
        paired = p.get("SEQUENCING")
        # Per sample paired, not implemented yet
        # pairedlist = p.get('SEQUENCING')
        # samplelist = p.get('SAMPLES')
        # x = samplelist.index(s.split(os.sep)[-1])
        # paired = pairedlist[x]
        ret.append(str(paired).replace(",", "_"))
    log.debug(logid + "PAIRED: " + str(ret))
    return str.join(",", ret)


@check_run
def checkstranded(sample: list, config: dict) -> str:
    """check strandedness of sample

    Parameters
    ----------
    sample : list
        samples
    config : dict
        config json file parsed by snakemake.common.configfile

    Returns
    -------
    str
        strandedness of sample
    """
    logid = scriptname + ".Params_checkstranded: "
    ret = list()
    stranded = ""
    for s in sample:
        check = conditiononly(s, config)
        log.debug(logid + "check: " + str(check))
        p = mu.sub_dict(config["SETTINGS"], check)
        log.debug(logid + "P: " + str(p))
        paired = p.get("SEQUENCING")
        # Per sample paired, not implemented yet
        # pairedlist = p.get('SEQUENCING')
        # samplelist = p.get('SAMPLES')
        # x = samplelist.index(s.split(os.sep)[-1])
        # paired = pairedlist[x]
        stranded = (
            paired.split(",")[1]
            if len(paired.split(",")) > 1 and paired.split(",")[1] != "unstranded"
            else ""
        )
    log.debug(logid + "STRANDEDNESS: " + str(stranded))
    return stranded


def resolve_sequencing(seqval: str, default_paired: str, default_stranded: str):
    """Resolve a per-workflow-step SEQUENCING override (e.g. 'single,rf') into
    a (paired, stranded) tuple, falling back to the given defaults (usually
    the global SETTINGS.SEQUENCING derived values, or the resolved MAPPING
    values for steps downstream of MAPPING) when the step itself does not
    define its own override. This allows e.g. a quantseq MAPPING step to run
    in 'single' mode while TRIMMING/DEDUP still see the global 'paired' mode
    needed for UMI extraction from read2.

    Parameters
    ----------
    seqval : str
        SEQUENCING value configured for the step, e.g. 'single,rf', or falsy
        if not configured for this step.
    default_paired : str
        paired value to fall back to when seqval is not set.
    default_stranded : str
        stranded value to fall back to when seqval is not set.

    Returns
    -------
    tuple(str, str)
        (paired, stranded)
    """
    if seqval and seqval != "":
        parts = str(seqval).split(",")
        return parts[0], (parts[1] if len(parts) > 1 else "")
    return default_paired, default_stranded


@check_run
def set_pairing(samples: list, config: dict) -> list:
    """get pairs of samples

    Parameters
    ----------
    samples : list
        samples
    config : dict
        config json file parsed by snakemake.common.configfile

    Returns
    -------
    list
        pairs of samples
    """
    logid = scriptname + ".Params_set_pairing: "
    ret = list()
    cond = conditiononly(samples[0], config)
    pconf = mu.sub_dict(config["PEAKS"], cond)
    log.debug(logid + "SAMPLES: " + str(samples))
    pairlist = pconf.get("COMPARABLE", config["PEAKS"].get("COMPARABLE"))
    log.debug(logid + "PAIRLIST: " + str(pairlist))
    if pairlist:
        for k, v in pairlist.items():
            for x in samples:
                if str(k) == str(os.path.basename(x)):
                    ret.extend(samplecond([x], config))
    else:
        return samples
    log.debug(logid + "return: " + str(ret))
    return ret


@check_run
def get_pairing(
    sample: str,
    stype: str,
    config: dict,
    samples: list,
    scombo: str = "",
    mode: str = "smk",
) -> str:
    """get pairs of samples

    Parameters
    ----------
    sample : str
        sample name
    stype : str
        unique or sorted
    config : dict
        config json file parsed by snakemake.common.configfile
    samples : list
        samples
    scombo : str, optional
        condition and tools to work on, by default ""
    mode : str, optional
        snakemake or nextflow mode, by default "smk"

    Returns
    -------
    str
        pairs of samples
    """
    logid = scriptname + ".Params_get_pairing: "
    cond = conditiononly(sample, config)
    pconf = mu.sub_dict(config["PEAKS"], cond)
    pairlist = pconf.get("COMPARABLE", config["PEAKS"].get("COMPARABLE"))
    matching = ""
    log.debug(
        logid
        + "PAIRLIST: "
        + str(pairlist)
        + " SAMPLE: "
        + str(sample)
        + " SAMPLES: "
        + str(samples)
        + " Combo: "
        + str(scombo)
    )
    if pairlist:
        for k, v in pairlist.items():
            if str(k) == str(os.path.basename(sample)):
                for x in samples:
                    if str(v) == str(os.path.basename(x)) and x != sample:
                        log.debug(logid + "Match found: " + str(v) + " : " + str(x))
                        try:
                            matching = samplecond([x], config)[0].replace("MAPPED/", "")
                        except (KeyError, IndexError):
                            matching = x
                        log.info(logid + "PAIRINGS: " + sample + ": " + str(matching))
        if not matching or matching == "":
            log.error(
                logid
                + f"COMPARABLE set in config but no fitting pair could be found for sample {sample} in {pairlist}. Please check config."
            )
        else:
            if mode == "smk":
                retstr = (
                    f"-c {os.sep.join(['MAPPED', scombo, matching])}_mapped_{stype}.bam"
                )
            else:
                retstr = f"{os.sep.join(['MAPPED', scombo, sample])}_mapped_{stype}.bam, {os.sep.join(['MAPPED', scombo, matching])}_mapped_{stype}.bam"

            log.debug(logid + retstr)
            return retstr
    else:
        log.warning(logid + "No matching sample found")
        return ""


@check_run
def post_checkpaired(sample: list, config: dict) -> str:
    """check if sample is paired

    Parameters
    ----------
    sample : list
        samples
    config : dict
        config json file parsed by snakemake.common.configfile

    Returns
    -------
    str
        pairedness of sample
    """
    logid = scriptname + ".Params_checkpaired: "
    ret = list()
    paired = ""
    for s in sample:
        log.debug(logid + "SAMPLE: " + str(sample))
        check = conditiononly(sample, config)
        p = mu.sub_dict(config["SETTINGS"], check)
        log.debug(logid + "P: " + str(p))
        paired = p.get("SEQUENCING").split(",")[0]
    log.debug(logid + "PAIRED: " + str(paired))
    return paired


@check_run
def check_IP(sample: list, config: dict) -> str:
    """get IP information from config

    Parameters
    ----------
    sample : list
        samples
    config : dict
        config json file parsed by snakemake.common.configfile

    Returns
    -------
    str
        IP information
    """
    logid = scriptname + ".Params_check_IP: "
    ret = list()
    ip = ""
    for s in sample:
        log.debug(logid + "SAMPLE: " + str(s))
        check = os.path.dirname(s).split(os.sep)
        tmplist = check
        for r in runstate_from_sample([s], config):
            if r not in tmplist:
                tmplist.extend(r)
        log.debug(logid + "TMP: " + str(tmplist))
        check = mu.get_from_dict(config["PEAKS"], tmplist)[0]
        log.debug(logid + "CHECK: " + str(check))
        if "IP" in check:
            ip = check["IP"]
        else:
            log.debug(logid + "Key IP not found in config")
    log.debug(logid + "IP is: " + str(ip))
    return str(ip)


@check_run
def check_tool_params(
    sample: str, runstate: str, config: dict, subconf: str, idx: str
) -> str:
    """check if tool parameters are set in config

    Parameters
    ----------
    sample : str
        sample name
    runstate : str
        condition tree part
    config : dict
        config json file parsed by snakemake.common.configfile
    subconf : str
        part of config
    idx : str
        tool step

    Returns
    -------
    str
        tool and sample specific parameters
    """
    try:
        par = tool_params(sample, runstate, config, subconf)["OPTIONS"][idx]
        if par != "":
            return par
        elif subconf == "MAPPING":
            return "std"
        else:
            return ""
    except (KeyError, IndexError, TypeError):
        if subconf == "MAPPING":
            return "std"
        else:
            return ""


@check_run
def comparable_as_string(config: dict, subwork: str) -> str:
    """get comparables as string from config

    Parameters
    ----------
    config : dict
        config json file parsed by snakemake.common.configfile
    subwork : str
        workflow part

    Returns
    -------
    str
        DE comparables
    """
    logid = scriptname + ".comparable_as_string: "
    log.debug(logid + "this is the config: " + str(config))
    check = config[subwork].get("COMPARABLE")
    if check:
        log.debug(logid + "determine comparables in " + subwork)
        complist = []
        compdict = config[subwork]["COMPARABLE"]
        for contrast in compdict:
            As = ""
            Bs = ""
            for condition in compdict[contrast][0]:
                As = (As + "+" + condition).strip("+")
            for condition in compdict[contrast][1]:
                Bs = (Bs + "+" + condition).strip("+")
            complist.append(f"{contrast}:{As}-vs-{Bs}")
        compstr = ",".join(complist)
        return compstr
    else:
        log.warning(
            logid + "no comparables found in " + subwork + ". Compare All vs. All."
        )
        groups_by_condition = list(mu.yield_from_dict("GROUPS", config))
        log.debug(logid + "Groups by condition: " + str(groups_by_condition))
        flattened = sorted(
            set(val for sublist in groups_by_condition for val in sublist)
        )
        log.debug(logid + "Flattened groups: " + str(flattened))
        combined = list(set(itertools.permutations(flattened, 2)))
        log.debug(logid + "Combined groups: " + str(combined))
        complist = []
        for key, value in combined:
            complist.append(f"{key}-vs-{value}:{key}-vs-{value}")
        compstr = ",".join(complist)
        log.debug(logid + "Comparables string: " + compstr)
        return compstr


@check_run
def get_combo_name(combinations: list) -> mu.NestedDefaultDict:
    """get combination of environments and tools

    Parameters
    ----------
    combinations : list
        list of conditions

    Returns
    -------
    mu.NestedDefaultDict
        dict of environments and tools per condition
    """
    logid = scriptname + ".Params_get_combo_name: "
    combname = mu.NestedDefaultDict()

    for condition in combinations:
        log.debug(logid + "CONDITION: " + str(condition))
        combname[condition]["envs"] = list()
        combname[condition]["works"] = list()
        combos = combinations[condition]

        for combi in combos:
            envs = list()
            works = list()
            for step in combi:
                if isinstance(step, list):
                    for substep in step:
                        for work, env in substep.items():
                            envs.append(env)
                            works.append(work)
                else:
                    for work, env in step.items():
                        envs.append(env)
                        works.append(work)
            combname[condition]["envs"].append(str.join("-", envs))
            combname[condition]["works"].append(str.join("-", works))

    log.debug(logid + "ComboName: " + str(combname))
    return combname


@check_run
def fixRunParameters(
    config: dict,
    env: str,
    sample: str,
    runstate: str,
    workflow: str,
    step: str,
    paramToFix: str,
    fixedParam: str,
) -> str:
    para = tool_params(sample, runstate, config, workflow, env)["OPTIONS"].get(step, "")
    if paramToFix in para and fixedParam:
        return re.sub(paramToFix + r" [0-9a-zA-Z\.\_\-\\/]+", fixedParam, para)
    else:
        log.warning(
            f"No parameter to fix {paramToFix} in INDEX options or fix {fixedParam} defined for {para}!"
        )
        return para


#
# Params.py ends here
