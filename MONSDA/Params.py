# Params.py ---
#
# Filename: Params.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Tue Sep 18 15:39:06 2018 (+0200)
# Version:
# Package-Requires: ()
# Last-Updated: Thu Feb  4 18:01:07 2021 (+0100)
#           By: Joerg Fallmann
#     Update #: 2888
# URL:
# Doc URL:
# Keywords:
# Compatibility:
#
#

# Commentary:
#
#
#

# Change Log:
#
#
#
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>.
#
#

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

import os
import sys
import re
import glob
import shutil
import json
import numpy as np
import heapq
import itertools
from operator import itemgetter
from natsort import natsorted
import traceback as tb
from io import StringIO
from Bio import SeqIO
import gzip
import inspect
import subprocess
import collections
from collections import defaultdict, OrderedDict
import six
import logging
import hashlib
from snakemake import load_configfile
import functools
import datetime


# cmd_subfolder = [os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"../MONSDA/lib"), os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"MONSDA/lib"), os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"../lib"), os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"lib")]
# for x in cmd_subfolder:
# if x not in sys.path:
#        sys.path.insert(0, x)

from MONSDA.Utils import *

try:
    scriptname = os.path.basename(inspect.stack()[-1].filename).replace(".py", "")
    log = logging.getLogger(scriptname)

    lvl = log.level if log.level else "INFO"
    for handler in log.handlers[:]:
        handler.close()
        log.removeHandler(handler)

    handler = logging.FileHandler("LOGS/MONSDA.log", mode="a")
    handler.setFormatter(
        logging.Formatter(
            fmt="%(asctime)s %(levelname)-8s %(name)-12s %(message)s",
            datefmt="%m-%d %H:%M",
        )
    )
    log.addHandler(handler)
    handler = logging.StreamHandler()
    handler.setFormatter(
        logging.Formatter(
            fmt="%(asctime)s %(levelname)-8s %(name)-12s %(message)s",
            datefmt="%m-%d %H:%M",
        )
    )
    log.addHandler(handler)
    log.setLevel(lvl)

except Exception:
    exc_type, exc_value, exc_tb = sys.exc_info()
    tbe = tb.TracebackException(
        exc_type,
        exc_value,
        exc_tb,
    )
    print("".join(tbe.format()), file=sys.stderr)


@check_run
def get_samples(config):
    logid = scriptname + ".Params_get_samples: "
    SAMPLES = [os.path.join(x) for x in sampleslong(config)]
    log.debug(logid + "SAMPLES_LONG: " + str(SAMPLES))
    check = [
        os.path.join("FASTQ", str(x).replace(".fastq.gz", "") + "*.fastq.gz")
        for x in SAMPLES
    ]
    RETSAMPLES = list()
    for i in range(len(check)):
        s = check[i]
        log.debug(logid + "SEARCHING: " + s)
        paired = checkpaired([SAMPLES[i]], config)
        log.debug(logid + "PAIRED: " + str(paired))
        f = glob.glob(s)
        log.debug(logid + "SAMPLECHECK: " + str(f))
        if f:
            f = list(set([str.join(os.sep, s.split(os.sep)[1:]) for s in f]))
            if "paired" not in paired:
                RETSAMPLES.extend(
                    list(
                        set(
                            [
                                os.path.join(
                                    os.path.dirname(s),
                                    re.sub(
                                        r"_r1.fastq.gz|_R1.fastq.gz|_r2.fastq.gz|_R2.fastq.gz|.fastq.gz",
                                        "",
                                        os.path.basename(s),
                                    ),
                                )
                                for s in f
                            ]
                        )
                    )
                )
                log.debug(logid + "PAIREDSAMPLES: " + str(f))
            else:
                RETSAMPLES.extend([x.replace(".fastq.gz", "") for x in f])
        else:
            log.debug(logid + "SAMPLECHECK: " + str(f))
    log.debug(logid + "SAMPLETEST: " + str(RETSAMPLES))
    if len(RETSAMPLES) < 1:
        log.error(logid + "No samples found, please check config file")
        sys.exit()

    log.debug(logid + "SAMPLES: " + str(RETSAMPLES))
    return RETSAMPLES


@check_run
def get_samples_postprocess(config, subwork):
    logid = scriptname + ".Params_get_samples_postprocess: "
    SAMPLES = [
        os.path.join(x)
        for x in sampleslong(config)
        if len(getFromDict(config[subwork], conditiononly(x, config))) > 0
        and (
            not config[subwork].get("EXCLUDE")
            or os.path.basename(x) not in config[subwork]["EXCLUDE"]
        )
    ]  # only those samples where postprocessing steps are defined for in config, should we make this a per condition check?
    log.debug(logid + "SAMPLES_LONG: " + str(SAMPLES))
    check = [
        os.path.join("FASTQ", str(x).replace(".fastq.gz", "") + "*.fastq.gz")
        for x in SAMPLES
    ]
    RETSAMPLES = list()
    for i in range(len(check)):
        s = check[i]
        paired = checkpaired([SAMPLES[i]], config)
        log.debug(logid + "PAIRED: " + str(paired))
        f = glob.glob(s)
        log.debug(logid + "SAMPLECHECK: " + str(f))
        if f:
            f = list(set([str.join(os.sep, s.split(os.sep)[1:]) for s in f]))
            if "paired" in paired:
                RETSAMPLES.extend(
                    list(
                        set(
                            [
                                os.path.join(
                                    os.path.dirname(s),
                                    re.sub(
                                        r"_r1.fastq.gz|_R1.fastq.gz|_r2.fastq.gz|_R2.fastq.gz|.fastq.gz",
                                        "",
                                        os.path.basename(s),
                                    ),
                                )
                                for s in f
                            ]
                        )
                    )
                )
                log.debug(logid + "PAIREDSAMPLES: " + str(f))
            else:
                RETSAMPLES.extend([x.replace(".fastq.gz", "") for x in f])
        else:
            log.debug(logid + "SAMPLECHECK: " + str(f))
    log.debug(logid + "SAMPLETEST: " + str(RETSAMPLES))
    if len(RETSAMPLES) < 1:
        log.error(
            logid
            + "No samples found for "
            + str(subwork)
            + ", please check config file"
        )
        sys.exit()

    log.debug(logid + "SAMPLES: " + str(RETSAMPLES))
    return RETSAMPLES


@check_run
def check_samples(config):
    logid = scriptname + ".Params_check_samples: "
    tocheck = [os.path.join(x) for x in sampleslong(config, nocheck=True)]
    log.debug(logid + "SEARCHING: " + str(tocheck))
    check = [
        os.path.join("FASTQ", str(x).replace(".fastq.gz", "") + "*.fastq.gz")
        for x in tocheck
    ]
    for i in range(len(check)):
        s = check[i]
        log.debug(logid + "SEARCHING: " + s)
        f = glob.glob(s)
        log.debug(logid + "SAMPLECHECK: " + str(f))
        if f:
            continue
        else:
            return False
    return True


@check_run
def download_samples(config):
    logid = scriptname + ".Params_download_samples: "
    SAMPLES = [os.path.join(x) for x in sampleslong(config, nocheck=True)]
    log.debug(logid + "DOWNLOAD_SAMPLES_LONG: " + str(SAMPLES))
    return SAMPLES


@check_run
def basecall_samples(config):
    logid = scriptname + ".Params_basecall_samples: "
    SAMPLES = [os.path.join(x) for x in sampleslong(config)]
    log.debug(logid + "SAMPLES_LONG: " + str(SAMPLES))
    check = [
        os.path.join("RAW", str(x).replace(".fast5", "") + "*.fast5") for x in SAMPLES
    ]
    RETSAMPLES = list()
    for i in range(len(check)):
        s = check[i]
        log.debug(logid + "SEARCHING: " + s)
        f = glob.glob(s)
        log.debug(logid + "SAMPLECHECK: " + str(f))
        if f:
            f = list(set([str.join(os.sep, s.split(os.sep)[1:]) for s in f]))
            RETSAMPLES.extend([x.replace(".fast5", "") for x in f])
    log.debug(logid + "SAMPLETEST: " + str(RETSAMPLES))
    if len(RETSAMPLES) < 1:
        log.error(logid + "No samples found, please check config file")
        sys.exit()

    log.debug(logid + "SAMPLES: " + str(RETSAMPLES))
    return RETSAMPLES


@check_run
def get_conditions(config):
    logid = scriptname + ".Params_conditions: "
    ret = list()
    for k in keysets_from_dict(config["SETTINGS"], "SAMPLES"):
        ret.append(k)
    log.debug(logid + str(ret))
    return list(set(ret))


@check_run
def get_samples_from_dir(search, config, nocheck=None):  # CHECK
    logid = scriptname + ".Params_get_samples_from_dir: "
    samples = [
        x.replace(" ", "")
        for x in list(set(getFromDict(config["SETTINGS"], search)[0]["SAMPLES"]))
    ]
    log.debug(logid + f"Samples: {samples}, Search: {search}, Check: {nocheck}")
    if nocheck is not None:
        samples = [os.sep.join([os.sep.join(search[0:]), s]) for s in samples]
        return list(set(samples))
    for x in range(
        len(search), len(search) - 2, -1
    ):  # For arbitrary depth of ics we append subdirectories until samples are found, maximum of one setting additional to file path is allowed
        pat = os.sep.join(["FASTQ", os.sep.join(search[0:x]), "*.fastq.gz"])
        log.debug(logid + "REGEX: " + str(pat) + "\t" + "SAMPLES: " + str(samples))
        check = natsorted(glob.glob(pat), key=lambda y: y.lower())
        log.debug(logid + "check: " + str(check))
        if len(check) > 0:
            ret = list()
            clean = list()
            for (
                c
            ) in (
                check
            ):  # If sample fits glob pattern but is not actually in the part of the config we are looking at this needs to be checked as checkpaired returns None otherwise, e.g. Condition1 has Sample abc_R1.fastq and Condition2 has Sample abcd_R1.fastq
                log.debug(logid + "check: " + str(c))
                x = c.split(os.sep)[-1]
                for s in samples:
                    log.debug(logid + "x: " + str(x))
                    log.debug(logid + "sample: " + str(s))
                    if s + "_R" in x or s + ".fastq.gz" == x:
                        log.debug(
                            logid
                            + "FOUND: "
                            + s
                            + "_R"
                            + " or "
                            + s
                            + ".fastq.gz"
                            + " in "
                            + x
                        )
                        clean.append(c)
                        break
            log.debug(logid + "checkclean: " + str(clean))
            paired = checkpaired(
                [os.sep.join([os.sep.join(search), clean[0].split(os.sep)[-1]])], config
            )
            if paired is not None and "paired" in paired:
                log.debug(
                    logid
                    + "SEARCHING: "
                    + str(
                        [
                            os.sep.join(
                                [
                                    os.sep.join(os.path.dirname(s).split(os.sep)[1:]),
                                    re.sub(
                                        r"_r1.fastq.gz|_R1.fastq.gz|_r2.fastq.gz|_R2.fastq.gz|.fastq.gz",
                                        "",
                                        os.path.basename(s),
                                    ),
                                ]
                            )
                            for s in clean
                        ]
                    )
                )
                ret.extend(
                    list(
                        set(
                            [
                                os.sep.join(
                                    [
                                        os.sep.join(
                                            os.path.dirname(s).split(os.sep)[1:]
                                        ),
                                        re.sub(
                                            r"_r1.fastq.gz|_R1.fastq.gz|_r2.fastq.gz|_R2.fastq.gz|.fastq.gz",
                                            "",
                                            os.path.basename(s),
                                        ),
                                    ]
                                )
                                for s in clean
                            ]
                        )
                    )
                )
                log.debug(logid + "FOUND: " + str(ret))
                renamelist = [
                    re.sub(r"_r\d", lambda pat: pat.group(1).upper(), s) for s in ret
                ]
                for i in range(len(renamelist)):
                    if renamelist[i] != ret[i]:
                        log.warning(
                            "SAMPLE NAMES CONTAIN LOWER CASE r1/r2 INSTEAD OF R1/R2 FOR PAIRED END SEQUENCING, THEY WILL BE RENAMED"
                        )
                        os.rename(ret[i], renamelist[i])
            else:
                log.debug(
                    logid
                    + "SEARCHING: "
                    + str(
                        [
                            os.sep.join(s.split(os.sep)[1:]).replace(".fastq.gz", "")
                            for s in clean
                        ]
                    )
                )
                ret.extend(
                    [
                        os.sep.join(s.split(os.sep)[1:]).replace(".fastq.gz", "")
                        for s in clean
                    ]
                )

            log.debug(logid + "RETURN: " + str(ret))
            return list(set(ret))
    log.error(logid + "NO SAMPLES FOUND")
    return list()


@check_run
def sampleslong(config, nocheck=None):
    logid = scriptname + ".Params_sampleslong: "
    log.debug(logid + f"Check: {nocheck}")
    tosearch = list()
    ret = list()
    for k in keysets_from_dict(config["SETTINGS"], "SAMPLES"):
        tosearch.append(list(k))
    log.debug(logid + "tosearch: " + str(tosearch))
    for search in tosearch:
        ret.extend(get_samples_from_dir(search, config, nocheck))
    log.debug(logid + str(ret))
    return ret


@check_run
def get_placeholder(config):
    ret = list()
    if "PH" in (config):
        for x in config["PH"]:
            ret.append(str(x))
    else:
        ret.append("_")
    return ret


@check_run
def get_cutoff_as_string(config, subwork, cf):
    logid = scriptname + ".get_cutoff: "
    cutoff = (
        str(config[subwork]["CUTOFFS"].get(cf))
        if config[subwork].get("CUTOFFS")
        else ".05"
        if cf == "pval"
        else "1.5"
    )
    log.info(logid + "CUTOFFS: " + str(cf) + ":" + cutoff)
    return cutoff


@check_run
def get_summary_dirs(config):
    logid = scriptname + ".get_summary_dirs: "
    ret = list()
    for work, tools in config["WORKFLOWS"].items():
        for tool in tools:
            ret.append(f"{work}/{tool.upper()}")
    log.debug(logid + str(ret))
    return ret


@check_run
def get_summary_files(config):
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
def create_skeleton(runner, skeleton=None):
    logid = scriptname + ".Params_create_skeleton: "
    if skeleton:
        for subdir in ["SubSnakes", "RAW", "FASTQ", "LOGS", "TMP", "JOB"]:
            makeoutdir(subdir)
            sys.exit(
                "Skeleton directories created, please add files and rerun without --skeleton option"
            )
    else:
        for subdir in [runner, "LOGS", "TMP", "JOBS"]:
            makeoutdir(subdir)
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


@check_run
def tool_params(sample, runstate, config, subconf, tool=None):
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
    mp = subDict(config[subconf], x)[tool] if tool else subDict(config[subconf], x)
    log.debug(logid + "DONE: " + str(mp))
    return mp


@check_run
def setting_per_sample(sample, runstate, config, setting, subconf=None):
    logid = scriptname + ".Params_setting_per_sample: "
    log.debug(logid + "Samples: " + str(sample))
    set = None
    x = sample.split(os.sep)[2:-1]
    if runstate is None:
        runstate = runstate_from_sample([sample], config)[0]
    if runstate not in x:
        x.append(runstate)
    subsetting = subDict(config["SETTINGS"], x).get(setting)

    if setting == "ANNOTATION":  # Special case is annotation
        subsetting = subsetting.get(
            "GTF", subsetting.get("GFF")
        )  # by default GTF format will be used

    if subconf:  # check specific setting for workflow part
        subset = (
            config[subconf].get(setting)
            if config[subconf].get(setting)
            else subDict(config[subconf], x).get(setting)
        )

    # Define which final setting is returned
    set = subset if subset else setting

    return set


@check_run
def get_reps(samples, config, analysis):
    logid = scriptname + ".Params_get_reps: "
    log.debug(logid + "Samples: " + str(samples))
    ret = defaultdict(list)
    for sample in samples:
        scond = (
            sample.split(os.sep)[4:-1]
            if len(sample.split(os.sep)) > 5
            else sample.split(os.sep)[3:-1]
        )
        log.debug(logid + "WORKING ON: " + str(sample) + " CONDITION: " + str(scond))
        partconf = subDict(config["SETTINGS"], scond)
        log.debug(logid + "CONF: " + str(partconf))

        Ex = config[analysis].get("EXCLUDE")
        if Ex and sample.split(os.sep)[-1] in Ex:
            continue
        ret["reps"].append(sample)
        wcfile = sample.split(os.sep)[-1].replace("_mapped_sorted_unique.counts.gz", "")
        idx = partconf["SAMPLES"].index(wcfile)
        scond.append(sample.split(os.sep)[-1])
        ret["pairs"].append(checkpaired_rep([str.join(os.sep, scond)], config))
        ret["conds"].append(partconf["GROUPS"][idx])
        if "BATCHES" in partconf and len(partconf["BATCHES"]) >= idx:
            ret["batches"].append(str(partconf["BATCHES"][idx]).replace(",", "_"))
        else:
            ret["batches"].append("1")
        if "TYPES" in partconf and len(partconf["TYPES"]) >= idx:
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
def get_diego_samples(samples, config, analysis):
    logid = scriptname + ".Params_get_diego_samples: "
    log.debug(logid + "Samples: " + str(samples))
    ret = defaultdict(list)
    for sample in samples:
        scond = sample.split(os.sep)[4:-1]
        log.debug(logid + "WORKING ON: " + str(sample) + " CONDITION: " + str(scond))
        partconf = subDict(config[analysis], scond)
        log.debug(logid + "CONF: " + str(partconf))
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
def get_diego_groups(samples, config, analysis):
    logid = scriptname + ".Params_get_diego_groups: "
    log.debug(logid + "Samples: " + str(samples))
    ret = defaultdict(list)
    for sample in samples:
        scond = sample.split(os.sep)[4:-1]
        log.debug(logid + "WORKING ON: " + str(sample) + " CONDITION: " + str(scond))
        partconf = subDict(config[analysis], scond)
        log.debug(logid + "CONF: " + str(partconf))
        wcfile = str.join("-", sample.split(os.sep)[-4:]).replace(
            "_mapped_sorted_unique.counts.gz", ""
        )
        checkfile = sample.split(os.sep)[-1].replace(
            "_mapped_sorted_unique.counts.gz", ""
        )
        idx = partconf["REPLICATES"].index(checkfile)
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
def env_bin_from_config(samples, config, subconf):
    logid = scriptname + ".Params_env_bin_from_config: "
    s = samples[0].split(os.sep)[:-1]
    mb, me = [None, None]
    for k in getFromDict(config[subconf], s):
        mb, me = k["BIN"], k["ENV"]
    return mb, me


@check_run
def env_bin_from_config2(samples, config, subconf):
    logid = scriptname + ".Params_env_bin_from_config2: "
    for s in samples:
        log.debug(logid + "S: " + str(s))
        log.debug(logid + "C: " + str(conditiononly(s, config)))
        check = conditiononly(s, config)

        for k in getFromDict(config[subconf], check):
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
def env_bin_from_config3(config, subconf):
    logid = scriptname + ".Params_env_bin_from_config3: "
    envkey = subconf + "ENV"
    binkey = subconf + "BIN"
    me = config[envkey]
    mb = config[binkey]
    log.debug(logid + str([str(mb), str(me)]))
    return mb, me


@check_run
def sample_from_path(path):
    ret = str(os.path.join(os.path.split(str(path))[-1]))
    return ret


@check_run
def runstate_from_sample(sample, config):
    logid = scriptname + ".Params_runstate_from_sample: "
    ret = list()
    for s in sample:
        if len(s.split(os.sep)) < 2:
            s = samplecond(s, config)
        if len(getFromDict(config["SETTINGS"], s.split(os.sep))) < 1:
            s = os.path.dirname(s)
            n = s.split(os.sep)[-1]
        log.debug(logid + "SAMPLE: " + s)
        try:
            c = getFromDict(config["SETTINGS"], s.split(os.sep))[0]
        except:
            c = None
        log.debug(logid + "SETTINGS: " + str(c))
        if dict_inst(c):
            if not c.get("SAMPLES"):
                for k, v in c.items():
                    log.debug(
                        logid + "k: " + str(k) + ", v: " + str(v) + " c: " + str(c)
                    )
                    if dict_inst(v) and v.get("SAMPLES"):
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
def samplecond(
    sample, config
):  # takes list of sample names (including .fastq.gz) and returns a list with their conditions as directory path without fastq.gz ( ["condition/of/sample", ... ])
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
                if sname in getFromDict(config["SETTINGS"], tmp)[0].get("SAMPLES"):
                    tmplist.append(r)
        log.debug(logid + "TMPLIST: " + str(tmplist))
        paired = checkpaired([s], config)
        if "paired" in paired:  # subDict(config['SETTINGS'], tmplist)['SEQUENCING']:
            s = re.sub(r"_[r|R|][1|2]\.", "", s)
        r = os.sep.join(tmplist)
        if r not in s:
            ret.append(os.sep.join([r, os.path.basename(s)]))
        else:
            ret.append(os.sep.join([os.path.dirname(s), os.path.basename(s)]))
    log.debug(logid + "RETURN: " + str(ret))
    return ret


@check_run
def conditiononly(sample, config):
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
            if len(getFromDict(config["SETTINGS"], tmp)) > 0 and sname in getFromDict(
                config["SETTINGS"], tmp
            )[0].get("SAMPLES"):
                ret.append(r)
    log.debug(logid + "ret: " + str(ret))
    return ret


@check_run
def checkpaired(sample, config):
    logid = scriptname + ".Params_checkpaired: "
    paired = ""
    for s in sample:
        log.debug(logid + "SAMPLE: " + str(s))
        check = conditiononly(s, config)
        log.debug(logid + "CHECK: " + str(check))
        p = subDict(config["SETTINGS"], check)
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
def checkpaired_rep(sample, config):
    logid = scriptname + ".Params_checkpaired_rep: "
    log.debug(logid + "SAMPLE: " + str(sample))
    ret = list()
    for s in sample:
        check = conditiononly(s, config)
        p = subDict(config["SETTINGS"], check)
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
def checkstranded(sample, config):
    logid = scriptname + ".Params_checkstranded: "
    ret = list()
    stranded = ""
    for s in sample:
        check = conditiononly(s, config)
        log.debug(logid + "check: " + str(check))
        p = subDict(config["SETTINGS"], check)
        log.debug(logid + "P: " + str(p))
        paired = p.get("SEQUENCING")
        # Per sample paired, not implemented yet
        # pairedlist = p.get('SEQUENCING')
        # samplelist = p.get('SAMPLES')
        # x = samplelist.index(s.split(os.sep)[-1])
        # paired = pairedlist[x]
        stranded = paired.split(",")[1] if len(paired.split(",")) > 1 else ""
    log.debug(logid + "STRANDEDNESS: " + str(stranded))
    return stranded


@check_run
def set_pairing(samples, config):
    logid = scriptname + ".Params_set_pairing: "
    ret = list()
    cond = conditiononly(samples[0], config)
    pconf = subDict(config["PEAKS"], cond)
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
def get_pairing(sample, stype, config, samples, scombo=""):
    logid = scriptname + ".Params_get_pairing: "
    cond = conditiononly(sample, config)
    pconf = subDict(config["PEAKS"], cond)
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
                        except:
                            matching = x
                        log.info(logid + "PAIRINGS: " + sample + ": " + str(matching))
        if not matching or matching == "":
            log.error(
                logid
                + f"COMPARABLE set in config but no fitting pair could be found for sample {sample} in {pairlist}. Please check config."
            )
        else:
            retstr = (
                "-c MAPPED"
                + os.sep
                + str(scombo)
                + os.sep
                + str(matching)
                + "_mapped_"
                + str(stype)
                + ".bam"
            )

            log.debug(logid + retstr)
            return retstr
    else:
        log.warning(logid + "No matching sample found")
        return ""


@check_run
def post_checkpaired(sample, config):
    logid = scriptname + ".Params_checkpaired: "
    ret = list()
    paired = ""
    for s in sample:
        log.debug(logid + "SAMPLE: " + str(sample))
        check = conditiononly(sample, config)
        p = subDict(config["SETTINGS"], check)
        log.debug(logid + "P: " + str(p))
        paired = p.get("SEQUENCING").split(",")[0]
    log.debug(logid + "PAIRED: " + str(paired))
    return paired


@check_run
def check_IP(sample, config):
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
        check = getFromDict(config["PEAKS"], tmplist)[0]
        log.debug(logid + "CHECK: " + str(check))
        if "IP" in check:
            ip = check["IP"]
        else:
            log.debug(logid + "Key IP not found in config")
    log.debug(logid + "IP is: " + str(ip))
    return str(ip)


@check_run
def check_tool_params(sample, runstate, config, subconf, idx):
    try:
        par = tool_params(sample, runstate, config, subconf)["OPTIONS"][idx]
        if par != "":
            return par
        elif subconf == "MAPPING":
            return "std"
        else:
            return ""
    except:
        if subconf == "MAPPING":
            return "std"
        else:
            return ""


@check_run
def comparable_as_string(config, subwork):
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
        groups_by_condition = list(yield_from_dict("GROUPS", config))
        flattened = sorted(
            set(val for sublist in groups_by_condition for val in sublist)
        )
        combined = list(set(itertools.combinations(flattened, 2)))
        complist = []
        for key, value in combined:
            complist.append(f"{key}vs{value}:{key}-vs-{value}")
        compstr = ",".join(complist)
        return compstr


@check_run
def get_combo_name(combinations):
    logid = scriptname + ".Params_get_combo_name: "
    combname = NestedDefaultDict()

    for condition in combinations:
        log.debug(logid + "CONDITION: " + str(condition))
        combname[condition]["envs"] = list()
        combname[condition]["works"] = list()
        combos = combinations[condition]

        for combi in combos:
            envs = list()
            works = list()
            for step in combi:
                for work, env in step.items():
                    envs.append(env)
                    works.append(work)
            combname[condition]["envs"].append(str.join("-", envs))
            combname[condition]["works"].append(str.join("-", works))

    log.debug(logid + "ComboName: " + str(combname))
    return combname


#
# Params.py ends here
