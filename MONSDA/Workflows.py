# Workflows.py ---
#
# Filename: Workflows.py
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
import traceback as tb
from io import StringIO
from Bio import SeqIO
import gzip
import inspect
import subprocess
import collections
from collections import defaultdict, OrderedDict
from pkg_resources import parse_version
import six
import logging
import hashlib
from snakemake import load_configfile
import functools
import datetime

from MONSDA.Utils import *
from MONSDA.Params import *

try:
    pythonversion = f"python{str(sys.version_info.major)}.{str(sys.version_info.minor)}"
    installpath = os.path.dirname(__file__).replace(
        os.sep.join(["lib", pythonversion, "site-packages", "MONSDA"]), "share"
    )
except:
    installpath = os.path.cwd()

workflowpath = os.path.join(installpath, "MONSDA", "workflows")
envpath = os.path.join(installpath, "MONSDA", "envs") + os.sep
binpath = os.path.join(installpath, "MONSDA", "scripts")
condapath = re.compile(r'conda:\s+"')
logfix = re.compile(r'loglevel="INFO"')

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


# Code:All subs from here on

###Check MONSDA version
@check_run
def ns_check_version(v, r):
    logid = scriptname + ".MONSDA_check_version: "

    if parse_version(str(v)) == parse_version(str(r)):
        return True
    else:
        log.debug(f"Mismatch in installed version {v} and required version {r}")
        return shutil.which("MONSDA")


##############################
########Snakemake Subs########
##############################


@check_run
def get_combo(wfs, config, conditions):
    logid = scriptname + ".Params_get_combo: "
    log.debug(logid + str(wfs) + str(conditions))
    combos = NestedDefaultDict()

    if wfs is None or len(wfs) < 1:
        return None

    for condition in conditions:
        ret = list()
        for subwork in wfs:
            listoftools, listofconfigs = create_subworkflow(
                config, subwork, [condition]
            )
            if listoftools is None:
                log.warning(
                    logid
                    + "No entry fits condition "
                    + str(condition)
                    + " for processing step "
                    + str(subwork)
                )
                return None

            tools = list()
            for k, v in [toolenvs for toolenvs in listoftools]:
                tools.append({subwork: k})
                log.debug(
                    logid
                    + "Preparing "
                    + str(subwork)
                    + " for condition "
                    + str(condition)
                    + " with Tool: "
                    + str(tools)
                )
            ret.append(tools)

        combos[condition] = itertools.product(*ret)

    # no debug log here or iterator will be consumed
    return combos


@check_run
def get_processes(config):
    logid = scriptname + ".Params_get_processes: "

    preprocess = subworkflows = postprocess = []

    # Define workflow stages
    pre = ["QC", "FETCH", "BASECALL"]
    sub = ["TRIMMING", "MAPPING", "DEDUP", "QC"]
    post = [
        "COUNTING",
        "TRACKS",
        "PEAKS",
        "DE",
        "DEU",
        "DAS",
        "DTU",
        "CIRCS",
        "ANNOTATE",
    ]

    wfs = [x.replace(" ", "") for x in config["WORKFLOWS"].split(",")]

    if "WORKFLOWS" in config:
        log.debug(
            logid
            + "CONFIG-WORKFLOWS: "
            + str(wfs)
            + "\t"
            + str(pre)
            + "\t"
            + str(sub)
            + "\t"
            + str(post)
        )
        subworkflows = [str(x) for x in wfs if x in sub]
        log.debug(logid + "Sub: " + str(subworkflows))
        if len(subworkflows) == 0 or subworkflows[0] == "":
            subworkflows = []
        preprocess = [x for x in wfs if x in pre]
        if len(preprocess) == 0 or preprocess[0] == "":
            preprocess = None
        log.debug(
            logid
            + "Intermediate-WORKFLOWS: "
            + str([preprocess, subworkflows, postprocess])
        )

        if (
            subworkflows
            and any(w in subworkflows for w in ["TRIMMING", "MAPPING", "DEDUP"])
            and preprocess
            and "QC" in preprocess
        ):
            preprocess.remove("QC")

        if (
            preprocess
            and "QC" in preprocess
            and not any(w in subworkflows for w in ["TRIMMING", "MAPPING", "DEDUP"])
        ):
            subworkflows.remove("QC")

        postprocess = [x for x in wfs if x in post]
        if len(postprocess) == 0 or postprocess[0] == "":
            postprocess = []
    else:
        log.error("NO WORKFLOWS DEFINED, NOTHING TO DO!")
        sys.exit()

    if preprocess:
        try:
            all([config[x] or x == "" for x in preprocess])
        except KeyError:
            log.warning(
                logid
                + "Not all required preprocessing steps have configuration in the config file"
            )

    if subworkflows:
        try:
            all([config[x] or x == "TRIMMING" or x == "" for x in subworkflows])
        except KeyError:
            log.warning(
                logid
                + "Not all required subworkflows have configuration in the config file"
            )

    if postprocess:
        try:
            all([config[x] or x == "" for x in postprocess])
        except KeyError:
            log.warning(
                logid
                + "Not all required postprocessing steps have configuration in the config file"
            )

    log.debug(logid + "WORKFLOWS: " + str([preprocess, subworkflows, postprocess]))

    return [preprocess, subworkflows, postprocess]


@check_run
def create_subworkflow(config, subwork, conditions, envs=None, stage=None):
    logid = scriptname + ".Workflows_create_subworkflow: "
    log.debug(
        logid
        + f"config:{config}, subwork:{subwork}, condition:{conditions}, stage:{stage}, envs:{envs}"
    )
    toollist = list()
    configs = list()

    for condition in conditions:
        try:
            env = str(subDict(config[subwork], condition)[stage + "ENV"])
        except:
            if "TOOLS" not in config[subwork]:
                log.error(
                    "No tool environment found for "
                    + subwork
                    + "! Either key ENV or TOOLS must be set for "
                    + str(condition)
                    + "!"
                )
            env = ""
        try:
            exe = str(subDict(config[subwork], condition)[stage + "BIN"])
        except:
            if "TOOLS" not in config[subwork]:
                log.error(
                    "No tool binary found for "
                    + subwork
                    + "! Either key BIN or TOOLS must be set for "
                    + str(condition)
                    + "!"
                )
            exe = ""

        tempconf = NestedDefaultDict()
        if env != "" and exe != "":
            toollist.append([env, exe])
            tempconf[subwork + "ENV"] = env
            tempconf[subwork + "BIN"] = exe
        try:
            tempconf["MAXTHREADS"] = config["MAXTHREADS"]
            BINS = config.get("BINS")
            if not BINS:
                BINS = binpath
            tempconf["BINS"] = BINS
        except KeyError:
            exc_type, exc_value, exc_tb = sys.exc_info()
            tbe = tb.TracebackException(
                exc_type,
                exc_value,
                exc_tb,
            )
            log.error("".join(tbe.format()))
        try:
            for key in ["SETTINGS", subwork]:
                if len(getFromDict(config[subwork], condition)) < 1:
                    if any(
                        [subwork == x for x in ["QC", "DEDUP", "TRIMMING", "MAPPING"]]
                    ):
                        log.error(
                            logid
                            + "Keys "
                            + str(condition)
                            + " not defined for "
                            + str(key)
                        )
                    else:
                        log.warning(
                            logid
                            + "Keys "
                            + str(condition)
                            + " not defined for "
                            + str(key)
                            + ", will be removed from SAMPLES for this analysis"
                        )
                        toollist.append([None, None])
                        configs.append(None)
                else:
                    if key == "SETTINGS":
                        tempconf[key] = subSetDict(config[key], condition)
                        if (
                            config.get("DEDUP") and "DEDUP" in config["WORKFLOWS"]
                        ) or config.get("RUNDEDUP"):
                            tempconf["RUNDEDUP"] = "enabled"
                        continue
                    if ("TOOLS" in config[key] and env == "" and exe == "") and len(
                        getFromDict(config[key], condition)
                    ) != 0:  # env and exe overrule TOOLS
                        for k, v in config[key][
                            "TOOLS"
                        ].items():  # each tool will be added to the config
                            if envs and k in envs and k != "None":
                                toollist.append([k, v])
                                tempconf[key]["TOOLS"][k] = config[key]["TOOLS"][k]
                                tc = list(condition)
                                tc.append(k)
                                tempconf[key].merge(subSetDict(config[key], tc))
                            elif not envs:
                                toollist.append([k, v])
                                tempconf[key] = subSetDict(config[key], condition)
                    else:
                        tempconf[key] = subSetDict(config[key], condition)

                    if stage and stage == "POST":
                        tempconf["SAMPLES"] = get_samples_postprocess(config, key)
                    else:
                        tempconf["SAMPLES"] = subDict(config["SETTINGS"], condition)[
                            "SAMPLES"
                        ]

                    if any(
                        [
                            subwork == x
                            for x in [
                                "PEAKS",
                                "DE",
                                "DEU",
                                "DAS",
                                "DTU",
                                "COUNTING",
                                "TRACKS",
                                "CIRCS",
                            ]
                        ]
                    ):
                        if not config[subwork].get("TOOLS"):
                            tempconf[subwork] = subSetDict(config[subwork], condition)
                        if config[subwork].get("CUTOFFS"):
                            tempconf[subwork]["CUTOFFS"] = config[subwork]["CUTOFFS"]
                        if subwork == "COUNTING":
                            tempconf["COUNTING"]["FEATURES"] = config["COUNTING"][
                                "FEATURES"
                            ]
                        if "COMPARABLE" in config[subwork]:
                            tempconf[subwork]["COMPARABLE"] = config[subwork][
                                "COMPARABLE"
                            ]

        except KeyError:
            exc_type, exc_value, exc_tb = sys.exc_info()
            tbe = tb.TracebackException(
                exc_type,
                exc_value,
                exc_tb,
            )
            log.error("".join(tbe.format()))

        configs.append(tempconf) if tempconf.get("SETTINGS", False) else configs.append(
            None
        )

    log.debug(logid + str([toollist, configs]))

    return toollist, configs


@check_run
def make_pre(
    subwork,
    config,
    samples,
    conditions,
    subdir,
    loglevel,
    state="",
    subname=None,
    combinations=None,
):
    logid = scriptname + ".Workflows_make_pre: "
    log.debug(logid + "WORK: " + str(subwork))

    subjobs = list()
    jobs = list()
    condapath = re.compile(r'conda:\s+"')
    logfix = re.compile(r'loglevel="INFO"')

    if combinations:
        combname = get_combo_name(combinations)
        for condition in combname:
            worklist = combname[condition]["works"]
            envlist = combname[condition]["envs"]
            subconf = NestedDefaultDict()
            add = list()

            smkf = os.path.abspath(os.path.join(workflowpath, "header.smk"))
            with open(smkf, "r") as smk:
                for line in smk.readlines():
                    line = re.sub(logfix, "loglevel='" + loglevel + "'", line)
                    line = re.sub(condapath, 'conda:  "' + envpath, line)
                    if "include: " in line:
                        line = fixinclude(
                            line, loglevel, condapath, envpath, workflowpath, logfix
                        )
                    add.append(line)

            for i in range(len(worklist)):
                log.debug(
                    logid + " SUBLISTS: " + str(worklist[i]) + "\t" + str(envlist[i])
                )
                works = worklist[i].split("-")
                envs = envlist[i].split("-")
                subjobs = list()

                # Add variable for combination string
                subjobs.append(
                    "\ncombo = '"
                    + str(envlist[i])
                    + "'\n\nwildcard_constraints:\n\tcombo = combo,\n \trawfile = '|'.join(list(SAMPLES)),\n\tfile = '|'.join(list(samplecond(SAMPLES, config))),\n\tread = \"R1|R2\"\n"
                )
                subjobs.append("\n\n")

                # Add rulethemall based on chosen workflows
                subjobs.append(
                    "".join(
                        rulethemall(
                            subwork,
                            config,
                            loglevel,
                            condapath,
                            logfix,
                            envlist[i],
                        )
                    )
                )

                for j in range(len(works)):
                    listoftools, listofconfigs = create_subworkflow(
                        config, works[j], [condition], envs
                    )

                    if listoftools is None:
                        log.warning(
                            logid
                            + "No entry fits condition "
                            + str(condition)
                            + " for processing step "
                            + str(works[j])
                        )
                        return None

                    sconf = listofconfigs[0]
                    sconf.pop("RUNDEDUP", None)  # cleanup
                    sconf.pop("PREDEDUP", None)  # cleanup

                    for a in range(0, len(listoftools)):
                        toolenv, toolbin = map(str, listoftools[a])
                        if toolenv != envs[j] or toolbin is None:
                            continue
                        sconf[works[j] + "ENV"] = toolenv
                        sconf[works[j] + "BIN"] = toolbin
                        subconf.update(sconf)
                        subname = toolenv + ".smk"

                        if works[j] == "QC":
                            subname = toolenv + "_raw.smk"

                        smkf = os.path.abspath(os.path.join(workflowpath, subname))
                        with open(smkf, "r") as smk:
                            for line in smk.readlines():
                                line = re.sub(
                                    logfix, "loglevel='" + loglevel + "'", line
                                )
                                line = re.sub(condapath, 'conda: "' + envpath, line)
                                if "include: " in line:
                                    line = fixinclude(
                                        line,
                                        loglevel,
                                        condapath,
                                        envpath,
                                        workflowpath,
                                        logfix,
                                    )
                                subjobs.append(line)
                            subjobs.append("\n\n")

                smkf = os.path.abspath(os.path.join(workflowpath, "footer.smk"))
                with open(smkf, "r") as smk:
                    for line in smk.readlines():
                        line = re.sub(logfix, "loglevel='" + loglevel + "'", line)
                        line = re.sub(condapath, 'conda: "' + envpath, line)
                        if "include: " in line:
                            line = fixinclude(
                                line, loglevel, condapath, envpath, workflowpath, logfix
                            )
                        subjobs.append(line)
                    subjobs.append("\n\n")

                smko = os.path.abspath(
                    os.path.join(
                        subdir,
                        "_".join(
                            [
                                "_".join(condition),
                                state + works[j],
                                toolenv,
                                "subsnake.smk",
                            ]
                        ),
                    )
                )
                if os.path.exists(smko):
                    os.rename(smko, smko + ".bak")
                with open(smko, "w") as smkout:
                    smkout.write("".join(add))
                    smkout.write("".join(subjobs))

                confo = os.path.abspath(
                    os.path.join(
                        subdir,
                        "_".join(
                            [
                                "_".join(condition),
                                state + subwork,
                                toolenv,
                                "subconfig.json",
                            ]
                        ),
                    )
                )
                if os.path.exists(confo):
                    os.rename(confo, confo + ".bak")
                with open(confo, "a") as confout:
                    json.dump(subconf, confout)

                jobs.append([smko, confo])

    else:
        for condition in conditions:
            subjobs = list()
            add = list()
            smkf = os.path.abspath(os.path.join(workflowpath, "header.smk"))
            with open(smkf, "r") as smk:
                for line in smk.readlines():
                    line = re.sub(logfix, "loglevel='" + loglevel + "'", line)
                    line = re.sub(condapath, 'conda: "' + envpath, line)
                    if "include: " in line:
                        line = fixinclude(
                            line, loglevel, condapath, envpath, workflowpath, logfix
                        )
                    add.append(line)
                add.append("\n\n")

            # Add variable for combination string
            add.append(
                'combo="'
                + os.sep
                + "\"\n\nwildcard_constraints:\n\tcombo = combo,\n \trawfile = '|'.join(list(SAMPLES)),\n\tfile = '|'.join(list(samplecond(SAMPLES, config))),\n\tread = \"R1|R2\"\n"
            )
            add.append("\n\n")

            listoftools, listofconfigs = create_subworkflow(
                config, subwork, [condition]
            )
            if listoftools is None:
                log.warning(
                    logid
                    + "No entry fits condition "
                    + str(condition)
                    + " for processing step "
                    + str(subwork)
                )
                return None

            sconf = listofconfigs[0]
            sconf.pop("RUNDEDUP", None)  # cleanup
            sconf.pop("PREDEDUP", None)  # cleanup

            if sconf is None:
                continue
            for i in range(0, len(listoftools)):
                toolenv, toolbin = map(str, listoftools[i])
                log.debug(logid + f"Env:{toolenv}, Bin:{toolbin}")
                if toolenv is None or toolbin is None:
                    continue
                subconf = NestedDefaultDict()
                sconf[subwork + "ENV"] = toolenv
                sconf[subwork + "BIN"] = toolbin
                subconf.update(sconf)
                subname = toolenv + ".smk"

                if subwork == "QC":
                    subname = toolenv + "_raw.smk"

                # Add rulethemall based on chosen workflows
                add.append(
                    "".join(
                        rulethemall(
                            subwork,
                            config,
                            loglevel,
                            condapath,
                            logfix,
                            subname.replace(".smk", ""),
                        )
                    )
                )  # RuleThemAll for snakemake depending on chosen workflows

                smkf = os.path.abspath(os.path.join(workflowpath, subname))
                with open(smkf, "r") as smk:
                    for line in smk.readlines():
                        line = re.sub(logfix, "loglevel='" + loglevel + "'", line)
                        line = re.sub(condapath, 'conda: "' + envpath, line)
                        if "include: " in line:
                            line = fixinclude(
                                line, loglevel, condapath, envpath, workflowpath, logfix
                            )
                        subjobs.append(line)
                    subjobs.append("\n\n")

                smkf = os.path.abspath(os.path.join(workflowpath, "footer.smk"))
                with open(smkf, "r") as smk:
                    for line in smk.readlines():
                        line = re.sub(logfix, "loglevel='" + loglevel + "'", line)
                        line = re.sub(condapath, 'conda: "' + envpath, line)
                        if "include: " in line:
                            line = fixinclude(
                                line, loglevel, condapath, envpath, workflowpath, logfix
                            )
                        subjobs.append(line)
                    subjobs.append("\n\n")

                smko = os.path.abspath(
                    os.path.join(
                        subdir,
                        "_".join(
                            [
                                "_".join(condition),
                                state + subwork,
                                toolenv,
                                "subsnake.smk",
                            ]
                        ),
                    )
                )

                if os.path.exists(smko):
                    os.rename(smko, smko + ".bak")
                with open(smko, "a") as smkout:
                    smkout.write(str.join("", add))
                    smkout.write("".join(subjobs))

                confo = os.path.abspath(
                    os.path.join(
                        subdir,
                        "_".join(
                            [
                                "_".join(condition),
                                state + subwork,
                                toolenv,
                                "subconfig.json",
                            ]
                        ),
                    )
                )
                if os.path.exists(confo):
                    os.rename(confo, confo + ".bak")
                with open(confo, "a") as confout:
                    json.dump(subconf, confout)

                jobs.append([smko, confo])

    log.debug(logid + "JOBS: " + str(jobs))
    return jobs


@check_run
def make_sub(
    subworkflows,
    config,
    samples,
    conditions,
    subdir,
    loglevel,
    subname=None,
    combinations=None,
):
    logid = scriptname + ".Workflows_make_sub: "

    log.info(logid + f"STARTING PROCESSING FOR {conditions}")

    jobs = list()
    condapath = re.compile(r'conda:\s+"')
    logfix = re.compile(r'loglevel="INFO"')

    if combinations:
        combname = get_combo_name(combinations)
        for condition in combname:
            worklist = combname[condition]["works"]
            envlist = combname[condition]["envs"]
            subconf = NestedDefaultDict()
            add = list()

            smkf = os.path.abspath(os.path.join(workflowpath, "header.smk"))
            with open(smkf, "r") as smk:
                for line in smk.readlines():
                    line = re.sub(logfix, "loglevel='" + loglevel + "'", line)
                    line = re.sub(condapath, 'conda: "' + envpath, line)
                    if "include: " in line:
                        line = fixinclude(
                            line, loglevel, condapath, envpath, workflowpath, logfix
                        )
                    add.append(line)

            for i in range(len(worklist)):
                log.debug(
                    logid + " SUBLISTS: " + str(worklist[i]) + "\t" + str(envlist[i])
                )
                works = worklist[i].split("-")
                envs = envlist[i].split("-")
                subjobs = list()

                # Add variable for combination string
                subjobs.append(
                    "\ncombo = '"
                    + str(envlist[i])
                    + "'\n\nwildcard_constraints:\n\tcombo = combo,\n \trawfile = '|'.join(list(SAMPLES)),\n\tfile = '|'.join(list(samplecond(SAMPLES, config))),\n\tread = \"R1|R2\"\n"
                )
                subjobs.append("\n\n")

                # Add rulethemall based on chosen workflows
                subjobs.append(
                    "".join(
                        rulethemall(
                            subworkflows,
                            config,
                            loglevel,
                            condapath,
                            logfix,
                            envlist[i],
                        )
                    )
                )

                for j in range(len(works)):
                    listoftools, listofconfigs = create_subworkflow(
                        config, works[j], [condition], envs
                    )

                    if listoftools is None:
                        log.warning(
                            logid
                            + "No entry fits condition "
                            + str(condition)
                            + " for processing step "
                            + str(works[j])
                        )
                        return None

                    sconf = listofconfigs[0]
                    for a in range(0, len(listoftools)):
                        toolenv, toolbin = map(str, listoftools[a])
                        if toolenv != envs[j] or toolbin is None:
                            continue
                        sconf[works[j] + "ENV"] = toolenv
                        sconf[works[j] + "BIN"] = toolbin

                        subconf.update(sconf)
                        subname = toolenv + ".smk"
                        log.debug(logid + f"SCONF:{sconf}, SUBCONF:{subconf}")
                        if (
                            works[j] == "QC"
                            and "TRIMMING" in works
                            and not "MAPPING" in works
                        ):
                            if "DEDUP" in works and "umitools" in envs:
                                subname = toolenv + "_dedup_trim.smk"
                            else:
                                subname = toolenv + "_trim.smk"

                        if (
                            works[j] == "QC"
                            and not "TRIMMING" in works
                            and not "MAPPING" in works
                        ):
                            if "DEDUP" in subworkflows and "umitools" in envs:
                                subname = toolenv + "_dedup.smk"
                            else:
                                subname = toolenv + "_raw.smk"

                        # Picard tools can be extended here
                        if works[j] == "DEDUP" and toolenv == "picard":
                            subname = toolenv + "_dedup.smk"
                            subconf.pop("PREDEDUP", None)
                        elif works[j] == "DEDUP" and toolenv == "umitools":
                            subconf["PREDEDUP"] = "enabled"

                        smkf = os.path.abspath(os.path.join(workflowpath, subname))
                        with open(smkf, "r") as smk:
                            for line in smk.readlines():
                                line = re.sub(condapath, 'conda: "' + envpath, line)
                                if "include: " in line:
                                    line = fixinclude(
                                        line,
                                        loglevel,
                                        condapath,
                                        envpath,
                                        workflowpath,
                                        logfix,
                                    )
                                subjobs.append(line)
                            subjobs.append("\n\n")

                if "MAPPING" in works:
                    smkf = os.path.abspath(os.path.join(workflowpath, "mapping.smk"))
                    with open(smkf, "r") as smk:
                        for line in smk.readlines():
                            line = re.sub(condapath, 'conda: "' + envpath, line)
                            if "include: " in line:
                                line = fixinclude(
                                    line,
                                    loglevel,
                                    condapath,
                                    envpath,
                                    workflowpath,
                                    logfix,
                                )
                            subjobs.append(line)
                        subjobs.append("\n\n")
                if "QC" in subworkflows:
                    smkf = os.path.abspath(os.path.join(workflowpath, "multiqc.smk"))
                    with open(smkf, "r") as smk:
                        for line in smk.readlines():
                            line = re.sub(condapath, 'conda: "' + envpath, line)
                            if "include: " in line:
                                line = fixinclude(
                                    line,
                                    loglevel,
                                    condapath,
                                    envpath,
                                    workflowpath,
                                    logfix,
                                )
                            subjobs.append(line)
                        subjobs.append("\n\n")

                # Append footer and write out subsnake and subconf per condition
                smkf = os.path.abspath(os.path.join(workflowpath, "footer.smk"))
                with open(smkf, "r") as smk:
                    for line in smk.readlines():
                        line = re.sub(condapath, 'conda: "' + envpath, line)
                        if "include: " in line:
                            line = fixinclude(
                                line, loglevel, condapath, envpath, workflowpath, logfix
                            )
                        subjobs.append(line)
                    subjobs.append("\n\n")

                smko = os.path.abspath(
                    os.path.join(
                        subdir,
                        "_".join(["_".join(condition), envlist[i], "subsnake.smk"]),
                    )
                )
                if os.path.exists(smko):
                    os.rename(smko, smko + ".bak")
                with open(smko, "w") as smkout:
                    smkout.write("".join(add))
                    smkout.write("".join(subjobs))

                confo = os.path.abspath(
                    os.path.join(
                        subdir,
                        "_".join(["_".join(condition), envlist[i], "subconfig.json"]),
                    )
                )
                if os.path.exists(confo):
                    os.rename(confo, confo + ".bak")
                with open(confo, "w") as confout:
                    json.dump(subconf, confout)

                jobs.append([smko, confo])

    else:
        for condition in conditions:
            subjobs = list()
            add = list()

            smkf = os.path.abspath(os.path.join(workflowpath, "header.smk"))
            with open(smkf, "r") as smk:
                for line in smk.readlines():
                    line = re.sub(logfix, "loglevel='" + loglevel + "'", line)
                    line = re.sub(condapath, 'conda: "' + envpath, line)
                    if "include: " in line:
                        line = fixinclude(
                            line, loglevel, condapath, envpath, workflowpath, logfix
                        )
                    add.append(line)
                add.append("\n\n")

            # Add variable for combination string
            add.append(
                "wildcard_constraints:\n\tcombo = combo,\n \trawfile = '|'.join(list(SAMPLES)),\n\tfile = '|'.join(list(samplecond(SAMPLES, config))),\n\tread = \"R1|R2\"\n"
            )
            add.append("\n\n")

            for subwork in subworkflows:
                log.debug(logid + "PREPARING " + str(subwork) + " " + str(condition))
                listoftools, listofconfigs = create_subworkflow(
                    config, subwork, [condition]
                )
                if listoftools is None:
                    log.warning(
                        logid
                        + "No entry fits condition "
                        + str(condition)
                        + " for processing step "
                        + str(subwork)
                    )
                    return None

                sconf = listofconfigs[0]
                # if any("umitools" in x for x in listoftools):
                #    sconf["PREDEDUP"] = "enabled"
                for i in range(0, len(listoftools)):
                    toolenv, toolbin = map(str, listoftools[i])
                    if toolenv is None or toolbin is None:
                        continue
                    subconf = NestedDefaultDict()
                    sconf[subwork + "ENV"] = toolenv
                    sconf[subwork + "BIN"] = toolbin
                    subconf.update(sconf)
                    subname = toolenv + ".smk"

                    if (
                        subwork == "QC"
                        and "TRIMMING" in subworkflows
                        and not "MAPPING" in subworkflows
                    ):
                        if "DEDUP" in subworkflows and not "picard" in any(
                            [x for x in listoftools]
                        ):
                            subname = toolenv + "_dedup_trim.smk"
                        else:
                            subname = toolenv + "_trim.smk"

                    # Picard tools can be extended here
                    if subwork == "DEDUP" and toolenv == "picard":
                        subname = toolenv + "_dedup.smk"
                        subconf.poop("PREDEDUP", None)
                    elif works[j] == "DEDUP" and toolenv == "umitools":
                        subconf["PREDEDUP"] = "enabled"
                    # Add rulethemall based on chosen workflows
                    add.append(
                        "".join(
                            rulethemall(
                                subwork,
                                config,
                                loglevel,
                                condapath,
                                logfix,
                                subname.replace(".smk", ""),
                            )
                        )
                    )

                    smkf = os.path.abspath(os.path.join(workflowpath, subname))
                    with open(smkf, "r") as smk:
                        for line in smk.readlines():
                            line = re.sub(condapath, 'conda: "' + envpath, line)
                            if "include: " in line:
                                line = fixinclude(
                                    line,
                                    loglevel,
                                    condapath,
                                    envpath,
                                    workflowpath,
                                    logfix,
                                )
                            subjobs.append(line)
                        subjobs.append("\n\n")

            if "MAPPING" in subworkflows:
                smkf = os.path.abspath(os.path.join(workflowpath, "mapping.smk"))
                with open(smkf, "r") as smk:
                    for line in smk.readlines():
                        line = re.sub(condapath, 'conda: "' + envpath, line)
                        if "include: " in line:
                            line = fixinclude(
                                line, loglevel, condapath, envpath, workflowpath, logfix
                            )
                        subjobs.append(line)
                    subjobs.append("\n\n")
                if "QC" in subworkflows:
                    smkf = os.path.abspath(os.path.join(workflowpath, "multiqc.smk"))
                    with open(smkf, "r") as smk:
                        for line in smk.readlines():
                            line = re.sub(condapath, 'conda: "' + envpath, line)
                            if "include: " in line:
                                line = fixinclude(
                                    line,
                                    loglevel,
                                    condapath,
                                    envpath,
                                    workflowpath,
                                    logfix,
                                )
                            subjobs.append(line)
                        subjobs.append("\n\n")

            smkf = os.path.abspath(os.path.join(workflowpath, "footer.smk"))
            with open(smkf, "r") as smk:
                for line in smk.readlines():
                    line = re.sub(logfix, "loglevel='" + loglevel + "'", line)
                    line = re.sub(condapath, 'conda: "' + envpath, line)
                    if "include: " in line:
                        line = fixinclude(
                            line, loglevel, condapath, envpath, workflowpath, logfix
                        )
                    subjobs.append(line)
                subjobs.append("\n\n")

            smko = os.path.abspath(
                os.path.join(subdir, "_".join(["_".join(condition), "subsnake.smk"]))
            )
            if os.path.exists(smko):
                os.rename(smko, smko + ".bak")
            with open(smko, "a") as smkout:
                smkout.write(str.join("", add))
                smkout.write(str.join("", subjobs))
            confo = os.path.abspath(
                os.path.join(subdir, "_".join(["_".join(condition), "subconfig.json"]))
            )
            if os.path.exists(confo):
                os.rename(confo, confo + ".bak")
            with open(confo, "a") as confout:
                json.dump(subconf, confout)

            jobs.append([smko, confo])

    return jobs


@check_run
def make_post(
    postworkflow,
    config,
    samples,
    conditions,
    subdir,
    loglevel,
    subname=None,
    combinations=None,
):
    logid = scriptname + ".Workflows_make_post: "

    log.debug(logid + f"STARTING POSTPROCESSING {postworkflow} FOR {conditions}")

    jobs = list()
    condapath = re.compile(r'conda:\s+"')
    logfix = re.compile(r'loglevel="INFO"')
    summary_tools_set = set()
    summary_tools_dict = dict()

    if combinations:
        combname = get_combo_name(combinations)
        subwork = postworkflow

        if subwork in ["DE", "DEU", "DAS", "DTU"]:

            condition = list(combname.keys())[0]
            envlist = combname[condition].get("envs")
            subconf = NestedDefaultDict()
            add = list()

            smkf = os.path.abspath(os.path.join(workflowpath, "header.smk"))
            with open(smkf, "r") as smk:
                for line in smk.readlines():
                    line = re.sub(logfix, "loglevel='" + loglevel + "'", line)
                    line = re.sub(condapath, 'conda: "' + envpath, line)
                    if "include: " in line:
                        line = fixinclude(
                            line, loglevel, condapath, envpath, workflowpath, logfix
                        )
                    add.append(line)

            for i in range(len(envlist)):
                envs = envlist[i].split("-")

                listoftools, listofconfigs = create_subworkflow(
                    config, subwork, combname, stage="POST"
                )

                if listoftools is None:
                    log.error(
                        logid + "No entry in config fits processing step" + str(subwork)
                    )

                sconf = listofconfigs[0]
                sconf.pop("PREDEDUP", None)  # cleanup

                for c in range(1, len(listofconfigs)):
                    sconf = merge_dicts(sconf, listofconfigs[c])

                for a in range(0, len(listoftools)):
                    subjobs = list()
                    toolenv, toolbin = map(str, listoftools[a])
                    for cond in combname.keys():
                        tc = list(cond)
                        tc.append(toolenv)
                        sconf[subwork] = merge_dicts(
                            sconf[subwork], subSetDict(config[subwork], tc)
                        )

                    if sconf[subwork].get("TOOLS"):
                        sconf[subwork]["TOOLS"] = subDict(
                            sconf[subwork]["TOOLS"], [toolenv]
                        )
                    if subwork in ["DE", "DEU", "DAS", "DTU"] and toolbin not in [
                        "deseq",
                        "diego",
                    ]:  # for all other postprocessing tools we have     more than     one         defined subworkflow
                        toolenv = toolenv + "_" + subwork

                    sconf[subwork + "ENV"] = toolenv
                    sconf[subwork + "BIN"] = toolbin

                    log.debug(
                        logid
                        + "POSTPROCESS: "
                        + str(subwork)
                        + " TOOL: "
                        + str(toolenv)
                    )

                    scombo = str(envlist[i]) if envlist[i] != "" else ""
                    combo = (
                        str.join(os.sep, [str(envlist[i]), toolenv])
                        if envlist[i] != ""
                        else toolenv
                    )

                    # Add variable for combination string
                    subjobs.append(
                        "\ncombo = '"
                        + combo
                        + "'\n"
                        + "\nscombo =      '"
                        + scombo
                        + "'\n"
                        + '\nwildcard_constraints:\n\tcombo = combo,\n\tscombo = scombo,\n\tread = "R1|R2",\n\ttype = "sorted|sorted_unique" if not rundedup else "sorted|sorted_unique|sorted_dedup|sorted_unique_dedup"'
                    )
                    subjobs.append("\n\n")
                    subconf.update(sconf)

                    subname = toolenv + ".smk"
                    smkf = os.path.abspath(os.path.join(workflowpath, subname))

                    with open(smkf, "r") as smk:
                        for line in smk.readlines():
                            line = re.sub(condapath, 'conda: "' + envpath, line)
                            if "include: " in line:
                                line = fixinclude(
                                    line,
                                    loglevel,
                                    condapath,
                                    envpath,
                                    workflowpath,
                                    logfix,
                                )
                            subjobs.append(line)
                        subjobs.append("\n\n")

                    # Append footer and write out subsnake and subconf per condition
                    smkf = os.path.abspath(os.path.join(workflowpath, "footer.smk"))
                    with open(smkf, "r") as smk:
                        for line in smk.readlines():
                            line = re.sub(condapath, 'conda: "' + envpath, line)
                            if "include: " in line:
                                line = fixinclude(
                                    line,
                                    loglevel,
                                    condapath,
                                    envpath,
                                    workflowpath,
                                    logfix,
                                )
                            subjobs.append(line)
                        subjobs.append("\n\n")

                    te = (
                        toolenv.split("_")[0] if "_" in toolenv else toolenv
                    )  # shorten toolenv if subwork is already added
                    smko = os.path.abspath(
                        os.path.join(
                            subdir,
                            "_".join(
                                [
                                    "_".join(condition),
                                    envlist[i],
                                    subwork,
                                    te,
                                    "subsnake.smk",
                                ]
                            ),
                        )
                    )
                    if os.path.exists(smko):
                        os.rename(smko, smko + ".bak")
                    with open(smko, "w") as smkout:
                        smkout.write("".join(add))
                        smkout.write("".join(subjobs))

                    confo = os.path.abspath(
                        os.path.join(
                            subdir,
                            "_".join(
                                [
                                    "_".join(condition),
                                    envlist[i],
                                    subwork,
                                    te,
                                    "subconfig.json",
                                ]
                            ),
                        )
                    )
                    if os.path.exists(confo):
                        os.rename(confo, confo + ".bak")
                    with open(confo, "w") as confout:
                        json.dump(subconf, confout)

                    jobs.append([smko, confo])

        else:
            for condition in combname:
                envlist = combname[condition].get("envs")
                log.debug(logid + f"POSTLISTS:{condition}, {subwork}, {envlist}")

                subconf = NestedDefaultDict()
                add = list()

                smkf = os.path.abspath(os.path.join(workflowpath, "header.smk"))
                with open(smkf, "r") as smk:
                    for line in smk.readlines():
                        line = re.sub(logfix, "loglevel='" + loglevel + "'", line)
                        line = re.sub(condapath, 'conda: "' + envpath, line)
                        if "include: " in line:
                            line = fixinclude(
                                line, loglevel, condapath, envpath, workflowpath, logfix
                            )
                        add.append(line)

                for i in range(len(envlist)):
                    envs = envlist[i].split("-")

                    listoftools, listofconfigs = create_subworkflow(
                        config, subwork, [condition], stage="POST"
                    )

                    if listoftools is None:
                        log.warning(
                            logid
                            + "No entry fits condition "
                            + str(condition)
                            + " for processing step "
                            + str(subwork)
                        )
                        continue

                    sconf = listofconfigs[0]
                    sconf.pop("PREDEDUP", None)  # cleanup

                    for c in range(1, len(listofconfigs)):
                        sconf = merge_dicts(sconf, listofconfigs[c])

                    for a in range(0, len(listoftools)):
                        subjobs = list()
                        toolenv, toolbin = map(str, listoftools[a])

                        if subwork == "CIRCS":
                            if toolenv == "ciri2" and "bwa" not in envs:
                                log.warning(
                                    "CIRI2 needs BWA mapped files, will skip input produced otherwise"
                                )
                                continue

                        tc = list(condition)
                        tc.append(toolenv)
                        sconf[subwork].update(subSetDict(config[subwork], tc))

                        if sconf[subwork].get("TOOLS"):
                            sconf[subwork]["TOOLS"] = subDict(
                                sconf[subwork]["TOOLS"], [toolenv]
                            )

                        sconf[subwork + "ENV"] = toolenv
                        sconf[subwork + "BIN"] = toolbin

                        log.debug(
                            logid
                            + "POSTPROCESS: "
                            + str(subwork)
                            + " CONDITION: "
                            + str(condition)
                            + " TOOL: "
                            + str(toolenv)
                        )

                        scombo = str(envlist[i]) if envlist[i] != "" else ""
                        combo = (
                            str.join(os.sep, [str(envlist[i]), toolenv])
                            if envlist[i] != ""
                            else toolenv
                        )

                        # Add variable for combination string
                        subjobs.append(
                            "\ncombo = '"
                            + combo
                            + "'\n"
                            + "\nscombo = '"
                            + scombo
                            + "'\n"
                            + '\nwildcard_constraints:\n\tcombo = combo,\n\tscombo = scombo,\n\tread = "R1|R2",\n\ttype = "sorted|sorted_unique" if not rundedup else "sorted|sorted_unique|sorted_dedup|sorted_unique_dedup"'
                        )
                        subjobs.append("\n\n")
                        subconf.update(sconf)

                        subname = toolenv + ".smk"
                        smkf = os.path.abspath(os.path.join(workflowpath, subname))

                        with open(smkf, "r") as smk:
                            for line in smk.readlines():
                                line = re.sub(condapath, 'conda: "' + envpath, line)
                                if "include: " in line:
                                    line = fixinclude(
                                        line,
                                        loglevel,
                                        condapath,
                                        envpath,
                                        workflowpath,
                                        logfix,
                                    )
                                subjobs.append(line)
                            subjobs.append("\n\n")

                        # Append footer and write out subsnake and subconf per condition
                        smkf = os.path.abspath(os.path.join(workflowpath, "footer.smk"))
                        with open(smkf, "r") as smk:
                            for line in smk.readlines():
                                line = re.sub(condapath, 'conda: "../', line)
                                if "include: " in line:
                                    line = fixinclude(
                                        line,
                                        loglevel,
                                        condapath,
                                        envpath,
                                        workflowpath,
                                        logfix,
                                    )
                                subjobs.append(line)
                            subjobs.append("\n\n")

                        te = (
                            toolenv.split("_")[0] if "_" in toolenv else toolenv
                        )  # shorten toolenv if subwork is already added
                        smko = os.path.abspath(
                            os.path.join(
                                subdir,
                                "_".join(
                                    [
                                        "_".join(condition),
                                        envlist[i],
                                        subwork,
                                        te,
                                        "subsnake.smk",
                                    ]
                                ),
                            )
                        )
                        if os.path.exists(smko):
                            os.rename(smko, smko + ".bak")
                        with open(smko, "w") as smkout:
                            smkout.write("".join(add))
                            smkout.write("".join(subjobs))

                        confo = os.path.abspath(
                            os.path.join(
                                subdir,
                                "_".join(
                                    [
                                        "_".join(condition),
                                        envlist[i],
                                        subwork,
                                        te,
                                        "subconfig.json",
                                    ]
                                ),
                            )
                        )
                        if os.path.exists(confo):
                            os.rename(confo, confo + ".bak")
                        with open(confo, "w") as confout:
                            json.dump(subconf, confout)

                        jobs.append([smko, confo])

    else:
        subwork = postworkflow
        add = list()

        for condition in conditions:
            subconf = NestedDefaultDict()
            smkf = os.path.abspath(os.path.join(workflowpath, "header.smk"))
            with open(smkf, "r") as smk:
                for line in smk.readlines():
                    line = re.sub(logfix, "loglevel='" + loglevel + "'", line)
                    line = re.sub(condapath, 'conda: "' + envpath, line)
                    if "include: " in line:
                        line = fixinclude(
                            line, loglevel, condapath, envpath, workflowpath, logfix
                        )
                    add.append(line)

            listoftools, listofconfigs = create_subworkflow(
                config, subwork, [condition], stage="POST"
            )

            if listoftools is None:
                log.warning(
                    logid
                    + "No entry fits condition "
                    + str(condition)
                    + " for processing step "
                    + str(subwork)
                )
                continue

            sconf = listofconfigs[0]
            sconf.pop("PREDEDUP", None)  # cleanup

            for a in range(0, len(listoftools)):
                subjobs = list()

                toolenv, toolbin = map(str, listoftools[a])
                if subwork in ["DE", "DEU", "DAS", "DTU"] and toolbin not in [
                    "deseq",
                    "diego",
                ]:  # for all other postprocessing tools we have more than one defined subworkflow
                    toolenv = toolenv + "_" + subwork
                    log.debug(logid + "toolenv: " + str(toolenv))

                sconf[subwork + "ENV"] = toolenv
                sconf[subwork + "BIN"] = toolbin

                scombo = ""
                combo = toolenv

                # Add variable for combination string
                subjobs.append(
                    "\ncombo = '"
                    + combo
                    + "'\n"
                    + "\nscombo = '"
                    + scombo
                    + "'\n"
                    + '\nwildcard_constraints:\n\tcombo = combo,\n\tscombo = scombo,\n\tread = "R1|R2",\n\ttype = "sorted|sorted_unique" if not rundedup else "sorted|sorted_unique|sorted_dedup|sorted_unique_dedup"'
                )
                subjobs.append("\n\n")
                subconf.update(sconf)

                subname = toolenv + ".smk"
                smkf = os.path.abspath(os.path.join(workflowpath, subname))

                with open(smkf, "r") as smk:
                    for line in smk.readlines():
                        line = re.sub(condapath, 'conda: "' + envpath, line)
                        if "include: " in line:
                            line = fixinclude(
                                line, loglevel, condapath, envpath, workflowpath, logfix
                            )
                        subjobs.append(line)
                    subjobs.append("\n\n")

                # Append footer and write out subsnake and subconf per condition
                smkf = os.path.abspath(os.path.join(workflowpath, "footer.smk"))
                with open(smkf, "r") as smk:
                    for line in smk.readlines():
                        line = re.sub(condapath, 'conda: "' + envpath, line)
                        if "include: " in line:
                            line = fixinclude(
                                line, loglevel, condapath, envpath, workflowpath, logfix
                            )
                        subjobs.append(line)
                    subjobs.append("\n\n")

                te = (
                    toolenv.split("_")[0] if "_" in toolenv else toolenv
                )  # shorten toolenv if subwork is already added
                smko = os.path.abspath(
                    os.path.join(
                        subdir,
                        "_".join(["_".join(condition), subwork, te, "subsnake.smk"]),
                    )
                )
                if os.path.exists(smko):
                    os.rename(smko, smko + ".bak")
                with open(smko, "w") as smkout:
                    smkout.write("".join(add))
                    smkout.write("".join(subjobs))

                confo = os.path.abspath(
                    os.path.join(
                        subdir,
                        "_".join(["_".join(condition), subwork, te, "subconfig.json"]),
                    )
                )
                if os.path.exists(confo):
                    os.rename(confo, confo + ".bak")
                with open(confo, "w") as confout:
                    json.dump(subconf, confout)

                jobs.append([smko, confo])

    return jobs


@check_run
def make_summary(config, subdir, loglevel, combinations=None):
    logid = scriptname + ".Workflows_make_summary: "

    output = "REPORTS/SUMMARY/summary.Rmd"
    jobs = list()
    lines = list()
    condapath = re.compile(r'conda:\s+"')
    logfix = re.compile(r'loglevel="INFO"')

    envlist = list()
    if combinations:
        combname = get_combo_name(combinations)
        for condition in combname:
            envlist = combname[condition]["envs"]

    # Add Header
    sum_path = os.path.join(installpath, "MONSDA", "scripts", "Analysis", "SUMMARY")
    rmd_header = os.path.abspath(os.path.join(sum_path, "header_summary.Rmd"))

    with open(rmd_header, "r") as read_file:
        for line in read_file.readlines():
            lines.append(line)
        lines.append("\n\n")

    # Add rMarkdown snippets
    if len(envlist) == 0:
        envlist.append("")
    for scombo in envlist:
        scombo = str(scombo) + os.sep
        snippath = str.join(os.sep, ["REPORTS", "SUMMARY", "RmdSnippets", scombo + "*"])
        snippets = [x for x in glob.glob(snippath) if os.path.isfile(x)]

        log.debug(f"{logid} snippets found: {snippets} for path: {snippath}")
        for snippet in snippets:
            with open(snippet, "r") as read_file:
                for line in read_file.readlines():
                    if line.startswith("# "):
                        lines = [[l] for l in "@$@".join(lines).split(line)]
                    lines[0].append(line)
            lines = "".join(sum(lines, [])).split("@$@")

    if os.path.exists(output):
        os.rename(output, output + ".bak")
    with open(output, "a") as writefile:
        for line in lines:
            writefile.write(line)
        writefile.write("\n\n")

    subjobs = list()

    smkf = os.path.abspath(os.path.join(workflowpath, "header.smk"))
    with open(smkf, "r") as smk:
        for line in smk.readlines():
            subjobs.append(line)
        subjobs.append("\n\n")

    smkf = os.path.abspath(os.path.join(workflowpath, "summary.smk"))
    with open(smkf, "r") as smk:
        for line in smk.readlines():
            line = re.sub(logfix, "loglevel='" + loglevel + "'", line)
            line = re.sub(condapath, 'conda: "' + envpath, line)
            if "include: " in line:
                line = fixinclude(
                    line, loglevel, condapath, envpath, workflowpath, logfix
                )
            subjobs.append(line)
        subjobs.append("\n\n")

    smkf = os.path.abspath(os.path.join(workflowpath, "footer.smk"))
    with open(smkf, "r") as smk:
        for line in smk.readlines():
            subjobs.append(line)
        subjobs.append("\n\n")

    smko = os.path.abspath(os.path.join(subdir, "summary_subsnake.smk"))
    if os.path.exists(smko):
        os.rename(smko, smko + ".bak")
    with open(smko, "a") as smkout:
        smkout.write("".join(subjobs))
        smkout.write("\n\n")

    subconf = NestedDefaultDict()
    for key in ["BINS", "MAXTHREADS", "SETTINGS"]:
        subconf[key] = config[key]

    confo = os.path.abspath(os.path.join(subdir, "summary_subconfig.json"))
    if os.path.exists(confo):
        os.rename(confo, confo + ".bak")
    with open(confo, "a") as confout:
        json.dump(subconf, confout)

    jobs.append([smko, confo])

    return jobs


@check_run
def rulethemall(subworkflows, config, loglevel, condapath, logfix, combo=""):
    logid = scriptname + ".Workflows_rulethemall: "

    allmap = (
        'rule themall:\n\tinput:\texpand("MAPPED/{combo}/{file}_mapped_{type}.bam", combo=combo, file=samplecond(SAMPLES, config), type=["sorted", "sorted_unique"])'
        if not "DEDUP" in subworkflows
        else 'rule themall:\n\tinput:\texpand("MAPPED/{combo}/{file}_mapped_{type}_dedup.bam", combo=combo, file=samplecond(SAMPLES, config), type=["sorted", "sorted_unique"])'
    )
    allqc = 'expand("QC/Multi/{combo}/{condition}/multiqc_report.html", condition=str.join(os.sep, conditiononly(SAMPLES[0], config)), combo=combo)'
    allrawqc = 'rule themall:\n\tinput:\texpand("QC/Multi{combo}{condition}/multiqc_report.html", condition=str.join(os.sep, conditiononly(SAMPLES[0], config)), combo=combo)'
    alltrimqc = 'rule themall:\n\tinput:\texpand("QC/Multi{combo}/{condition}/multiqc_report.html", condition=str.join(os.sep, conditiononly(SAMPLES[0], config)), combo=combo)'
    alltrim = 'rule themall:\n\tinput: expand("TRIMMED_FASTQ/{combo}/{file}_{read}_trimmed.fastq.gz", combo=combo, file=samplecond(SAMPLES, config), read=["R1","R2"]) if paired == \'paired\' else expand("TRIMMED_FASTQ/{combo}/{file}_trimmed.fastq.gz", combo=combo, file=samplecond(SAMPLES, config))'
    alldedupqc = 'rule themall:\n\tinput:\texpand("QC/Multi{combo}/{condition}/multiqc_report.html", combo=combo, condition=str.join(os.sep, conditiononly(SAMPLES[0], config)))'
    alldedup = 'rule themall:\n\tinput: expand("DEDUP_FASTQ/{combo}/{file}_{read}_dedup.fastq.gz", combo=combo, file=samplecond(SAMPLES, config), read=["R1","R2"]) if paired == \'paired\' else expand("DEDUP_FASTQ/{combo}/{file}_dedup.fastq.gz", combo=combo, file=samplecond(SAMPLES, config))'
    alltrimdedupqc = 'rule themall:\n\tinput:\texpand("QC/Multi/{combo}/{condition}/multiqc_report.html", condition=str.join(os.sep, conditiononly(SAMPLES[0], config)), combo=combo)'

    todos = list()

    if "QC" in subworkflows and "QC" in config:
        makeoutdir("QC")
        if "MAPPING" in subworkflows:
            todos.append(allmap + ",\n\t\t" + allqc + "\n\n")
        else:
            if "TRIMMING" in subworkflows and "DEDUP" not in subworkflows:
                todos.append(alltrimqc + "\n\n")
            elif "TRIMMING" in subworkflows and "DEDUP" in subworkflows:
                todos.append(alltrimdedupqc + "\n\n")
            elif "DEDUP" in subworkflows and "TRIMMING" not in subworkflows:
                todos.append(alldedupqc + "\n\n")
            else:
                todos.append(allrawqc + "\n\n")

    if "MAPPING" in subworkflows and "QC" not in subworkflows:
        log.debug(logid + "Mapping without QC!")
        todos.append(allmap + "\n\n")

    if "MAPPING" in subworkflows and "TRIMMING" not in subworkflows:
        log.debug(logid + "Simulated read trimming only!")
        makeoutdir("TRIMMED_FASTQ")
        smkf = os.path.abspath(os.path.join(workflowpath, "simulatetrim.smk"))
        with open(smkf, "r") as smk:
            for line in smk.readlines():
                line = re.sub(logfix, "loglevel='" + loglevel + "'", line)
                line = re.sub(condapath, 'conda:  "' + envpath, line)
                if "include: " in line:
                    line = fixinclude(
                        line, loglevel, condapath, envpath, workflowpath, logfix
                    )
                todos.append(line)
        todos.append("\n\n")

    if (
        "TRIMMING" in subworkflows
        and "QC" not in subworkflows
        and "MAPPING" not in subworkflows
    ):
        log.debug(logid + "Trimming without QC!")
        todos.append(alltrim + "\n\n")

    if (
        "DEDUP" in subworkflows
        and "QC" not in subworkflows
        and "TRIMMING" not in subworkflows
        and "MAPPING" not in subworkflows
    ):
        log.debug(logid + "DEDUP without QC!")
        todos.append(alldedup + "\n\n")

    return todos


@check_run
def fixinclude(
    line,
    loglevel,
    condapath=condapath,
    envpath=envpath,
    workflowpath=workflowpath,
    logfix=logfix,
):

    logid = scriptname + ".Workflows_fixinclude: "

    linelist = list()
    toinclude = str.split(line)[-1].replace('"', "")
    toinclude = str.join(os.sep, [workflowpath, toinclude])
    with open(toinclude, "r") as incl:
        for line in incl.readlines():
            line = re.sub(logfix, "loglevel='" + loglevel + "'", line)
            line = re.sub(condapath, 'conda:  "' + envpath, line)
            if "include: " in line:
                line = fixinclude(
                    line, loglevel, condapath, envpath, workflowpath, logfix
                )
            linelist.append(line)
    toinclude = str.join("", linelist)
    return toinclude


##############################
########Nextflow Subs########
##############################


@check_run
def nf_check_version(v):
    logid = scriptname + ".nf_check_version: "

    if shutil.which("nextflow"):
        jobtorun = ["nextflow", "-v"]
        out = subprocess.run(jobtorun, stdout=subprocess.PIPE)
        check = out.stdout.decode("utf-8").split(" ")[-1]
        if parse_version(v) < parse_version(check):
            log.debug(logid + check)
            return True
    else:
        return shutil.which("nextflow")


@check_run
def nf_fetch_params(
    configfile, condition=None, combi=None
):  # replaces header.smk for nextflow workflows
    logid = scriptname + ".nf_fetch_params: "
    log.debug(logid + str(configfile) + " " + str(condition) + " " + str(combi))

    config = load_configfile(configfile)
    retconf = collections.defaultdict()

    BINS = config["BINS"]
    if not BINS:
        BINS = binpath

    MAXTHREAD = int(config["MAXTHREADS"])
    SAMPLES = (
        [os.path.join(x) for x in sampleslong(config)]
        if not config.get("FETCH", False)
        else [os.path.join(x) for x in download_samples(config)]
        if not config.get("BASECALL", False)
        else [os.path.join(x) for x in basecall_samples(config)]
    )
    if len(SAMPLES) < 1:
        log.error(logid + "No samples found, please check config file")
        sys.exit(logid + "ERROR: No samples found, please check config file")

    SETUP = keysets_from_dict(config["SETTINGS"], "SAMPLES")[0]
    SETS = os.sep.join(SETUP)
    SETTINGS = subDict(config["SETTINGS"], SETUP)

    # Parse SETTINGS
    REFERENCE = SETTINGS.get("REFERENCE")
    REFDIR = str(os.path.dirname(REFERENCE))
    INDEX = SETTINGS.get("INDEX")
    INDEX2 = SETTINGS.get("INDEX")
    UIDX = SETTINGS.get("UIDX")
    PREFIX = SETTINGS.get("PREFIX")
    ANNO = SETTINGS.get("ANNOTATION")
    ANNOTATION = ANNO.get("GTF") if ANNO else ""
    IP = SETTINGS.get("IP")

    rundedup = config.get("RUNDEDUP")
    prededup = config.get("PREDEDUP")

    if rundedup:
        if prededup:
            log.debug("(PRE)DEDUPLICATION ENABLED")
        else:
            log.debug("DEDUPLICATION ENABLED")

    paired = checkpaired([SAMPLES[0]], config)
    if paired == "paired":
        log.debug("RUNNING NEXTFLOW IN PAIRED READ MODE")

    stranded = checkstranded([SAMPLES[0]], config)
    if stranded != "":
        log.debug("RUNNING NEXTFLOW WITH STRANDEDNESS " + str(stranded))

    # save in return dict
    retconf["BINS"] = BINS
    retconf["MAXTHREAD"] = MAXTHREAD
    retconf["SAMPLES"] = str.join(",", SAMPLES)
    LONGSAMPLES = samplecond(SAMPLES, config)
    retconf["LONGSAMPLES"] = str.join(",", LONGSAMPLES)
    SHORTSAMPLES = [os.path.basename(x) for x in SAMPLES]
    retconf["SHORTSAMPLES"] = str.join(",", SHORTSAMPLES)
    retconf["CONDITION"] = os.sep.join(condition) if condition else SETS
    if combi:
        retconf["COMBO"] = combi[0] + os.sep if combi[0] != "" else None
        retconf["SCOMBO"] = combi[1] if combi[1] else None

    sample = SAMPLES[0]
    lsample = LONGSAMPLES[0]
    if paired:
        retconf["PAIRED"] = paired
    if stranded:
        retconf["STRANDED"] = stranded
    if rundedup:
        retconf["RUNDEDUP"] = rundedup
    if prededup:
        retconf["PREDEDUP"] = prededup

    # MAPPING Variables
    if "MAPPING" in config:
        MAPCONF = subDict(config["MAPPING"], SETUP)
        MAPPERBIN, MAPPERENV = env_bin_from_config3(config, "MAPPING")
        MAPOPT = MAPCONF.get(MAPPERENV).get("OPTIONS")
        log.debug(logid + "MAPPINGCONFIG: " + str(SETUP) + "\t" + str(MAPCONF))
        REF = MAPCONF.get("REFERENCE", MAPCONF[MAPPERENV].get("REFERENCE"))
        MANNO = MAPCONF.get("ANNOTATION", MAPCONF[MAPPERENV].get("ANNOTATION"))
        if REF:
            REFERENCE = REF
            REFDIR = str(os.path.dirname(REFERENCE))
        if MANNO:
            ANNOTATION = MANNO
        else:
            ANNOTATION = (
                ANNO.get("GTF")
                if "GTF" in ANNO and ANNO.get("GTF") != ""
                else ANNO.get("GFF")
            )  # by default GTF format will be used
        PRE = MAPCONF.get(
            "PREFIX",
            MAPCONF.get("EXTENSION", MAPOPT.get("PREFIX", MAPOPT.get("EXTENSION"))),
        )
        if PRE and PRE is not None:
            PREFIX = PRE
        if not PREFIX or PREFIX is None:
            PREFIX = MAPPERENV
        IDX = MAPCONF.get("INDEX", MAPCONF[MAPPERENV].get("INDEX"))
        if IDX:
            INDEX = IDX
        if not INDEX:
            INDEX = str.join(os.sep, [REFDIR, "INDICES", MAPPERENV]) + ".idx"
        unikey = get_dict_hash(
            subDict(
                tool_params(SAMPLES[0], None, config, "MAPPING", MAPPERENV)["OPTIONS"],
                ["INDEX"],
            )
        )
        UIDX = f"{REFDIR}/INDICES/{MAPPERENV}_{unikey}"
        UIDXNAME = f"{MAPPERENV}_{unikey}"
        INDICES = INDEX.split(",") if INDEX else list(UIDX)
        INDEX = (
            str(os.path.abspath(INDICES[0]))
            if str(os.path.abspath(INDICES[0])) not in UIDX
            else str(os.path.abspath(INDICES[0])) + "_idx"
        )
        if len(INDICES) > 1:
            if str(os.path.abspath(INDICES[1])) not in UIDX:
                INDEX2 = str(os.path.abspath(INDICES[1]))
            else:
                INDEX2 = str(os.path.abspath(INDICES[1])) + "_idx"
        else:
            INDEX2 = None

        retconf["MAPPINGREF"] = REFERENCE
        retconf["MAPPINGREFDIR"] = REFDIR
        retconf["MAPPINGANNO"] = ANNOTATION
        retconf["MAPPINGIDX"] = INDEX
        retconf["MAPPINGIDX2"] = INDEX2
        retconf["MAPPINGUIDX"] = UIDX
        retconf["MAPPINGUIDXNAME"] = UIDXNAME
        retconf["MAPPINGPREFIX"] = PREFIX

    # Peak Calling Variables
    if "PEAKS" in config:
        PEAKCONF = subDict(config["PEAKS"], SETUP)
        REF = PEAKCONF.get("REFERENCE")
        if REF:
            REFERENCE = REF
            REFDIR = str(os.path.dirname(REFERENCE))
        ANNOPEAK = PEAKCONF.get("ANNOTATION")
        if ANNOPEAK:
            ANNOTATION = ANNOPEAK
        else:
            ANNOTATION = (
                ANNO.get("GTF") if "GTF" in ANNO else ANNO.get("GFF")
            )  # by default GTF forma
        if not IP:
            IP = check_IP(SAMPLES, config)
        log.debug(logid + "Running Peak finding for " + IP + " protocol")

        retconf["PEAKREF"] = REFERENCE
        retconf["PEAKREFDIR"] = REFDIR
        retconf["PEAKANNO"] = ANNOTATION
        retconf["PEAKIP"] = IP

    # TRACKS/COUNTING Variables
    for x in ["TRACKS", "COUNTING"]:
        if x in config:
            XCONF = subDict(config[x], SETUP)
            log.debug(logid + "XCONFIG: " + str(SETUP) + "\t" + str(XCONF))
            REF = XCONF.get("REFERENCE")
            XANNO = XCONF.get("ANNOTATION")
            if XANNO:
                ANNOTATION = XANNO
            else:
                ANNOTATION = (
                    ANNO.get("GTF") if "GTF" in ANNO else ANNO.get("GFF")
                )  # by default GTF forma
            if REF:
                REFERENCE = REF
                REFDIR = str(os.path.dirname(REFERENCE))
        retconf["TRACKSREF"] = REFERENCE
        retconf["TRACKSREFDIR"] = REFDIR
        retconf["TRACKSANNO"] = ANNOTATION

    # DE/DEU/DAS/DTU Variables
    for x in ["DE", "DEU", "DAS", "DTU"]:
        if x in config:
            XCONF = subDict(config[x], SETUP)
            log.debug(logid + "XCONFIG: " + str(SETUP) + "\t" + str(XCONF))
            REF = XCONF.get("REFERENCE")
            XANNO = XCONF.get("ANNOTATION")
            if XANNO:
                ANNOTATION = XANNO
            else:
                ANNOTATION = (
                    ANNO.get("GTF") if "GTF" in ANNO else ANNO.get("GFF")
                )  # by default GTF forma
            if REF:
                REFERENCE = REF
                REFDIR = str(os.path.dirname(REFERENCE))
        retconf["DEREF"] = REFERENCE
        retconf["DEREFDIR"] = REFDIR
        retconf["DEANNO"] = ANNOTATION

    # CIRCS Variables
    if "CIRCS" in config:
        CIRCCONF = subDict(config["CIRCS"], SETUP)
        log.debug(logid + "CIRCCONFIG: " + str(SETUP) + "\t" + str(CIRCCONF))
        REF = CIRCCONF.get("REFERENCE")
        if REF:
            REFERENCE = REF
            REFDIR = str(os.path.dirname(REFERENCE))
        CANNO = CIRCCONF.get("ANNOTATION")
        if CANNO:
            ANNOTATION = CANNO
        else:
            ANNOTATION = (
                ANNO.get("GTF") if "GTF" in ANNO else ANNO.get("GFF")
            )  # by default GTF format will be used
        retconf["CIRCREF"] = REFERENCE
        retconf["CIRCREFDIR"] = REFDIR
        retconf["CIRCANNO"] = ANNOTATION

    retconf["REFERENCE"] = REFERENCE
    retconf["REFDIR"] = REFDIR
    retconf["IP"] = IP

    return retconf


@check_run
def nf_tool_params(
    sample, runstate, config, subwork, toolenv, toolbin, workflows=None, condition=None
):
    logid = scriptname + ".nf_tool_params: "
    log.debug(
        logid
        + "Samples: "
        + str(sample)
        + ", "
        + str.join(
            ",",
            map(
                str, [runstate, config, subwork, toolenv, toolbin, workflows, condition]
            ),
        )
    )

    if " " in toolbin:
        toolbin = toolbin.replace(" ", "_")

    mp = OrderedDict()
    x = sample.split(os.sep)[:-1]
    if runstate is None:
        runstate = runstate_from_sample([sample], config)[0]
    if runstate not in x:
        x.append(runstate)

    tp = list()

    if not workflows:
        mp = (
            subDict(config[subwork], x)[toolenv]["OPTIONS"]
            if toolenv
            else subDict(config[subwork], x)["OPTIONS"]
        )
        tp.append(
            "--" + subwork + "ENV " + toolenv + " --" + subwork + "BIN " + toolbin + " "
        )

        toolpar = list()
        for key, val in mp.items():
            pars = val if val and val != "" else None
            if pars:
                tp.append("--" + toolenv + "_params_" + str(key) + " '" + pars + "'")
    else:
        for subwork in workflows:
            sd = subDict(config[subwork], condition)
            log.debug(logid + "SD: " + str(sd))
            if sd.get("ENV") and sd.get("BIN"):
                toolenv, toolbin = map(str, [sd["ENV"], sd["BIN"]])

            mp = sd[toolenv]["OPTIONS"] if toolenv else sd["OPTIONS"]
            tp.append(
                "--"
                + subwork
                + "ENV "
                + toolenv
                + " --"
                + subwork
                + "BIN "
                + toolbin
                + " "
            )

            toolpar = list()
            for key, val in mp.items():
                pars = val if val and val != "" else None
                if pars:
                    tp.append(
                        "--" + toolenv + "_params_" + str(key) + " '" + pars + "'"
                    )

    log.debug(logid + "DONE: " + str(tp))
    return " ".join(tp)


@check_run
def nf_get_processes(config):
    logid = scriptname + ".Workflows_nf_get_processes: "

    preprocess = subworkflows = postprocess = []

    # Define workflow stages
    pre = ["QC", "FETCH"]  # , 'BASECALL']
    sub = ["TRIMMING", "MAPPING", "QC", "DEDUP"]
    post = (
        []
    )  # ['COUNTING', 'TRACKS', 'PEAKS', 'DE', 'DEU', 'DAS', 'DTU', 'ANNOTATE']  # Not implemented yet, TODO

    wfs = [x.replace(" ", "") for x in config["WORKFLOWS"].split(",")]

    if "WORKFLOWS" in config:
        log.debug(
            logid
            + "CONFIG-WORKFLOWS: "
            + str(wfs)
            + "\t"
            + str(pre)
            + "\t"
            + str(sub)
            + "\t"
            + str(post)
        )
        subworkflows = [str(x) for x in wfs if x in sub]
        log.debug(logid + "Sub: " + str(subworkflows))
        if len(subworkflows) == 0 or subworkflows[0] == "":
            subworkflows = []
        preprocess = [x for x in wfs if x in pre]
        if len(preprocess) == 0 or preprocess[0] == "":
            preprocess = None
        log.debug(
            logid
            + "Intermediate-WORKFLOWS: "
            + str([preprocess, subworkflows, postprocess])
        )

        if (
            subworkflows
            and any(w in subworkflows for w in ["TRIMMING", "MAPPING", "DEDUP"])
            and preprocess
            and "QC" in preprocess
        ):
            preprocess.remove("QC")

        if (
            preprocess
            and "QC" in preprocess
            and not any(w in subworkflows for w in ["TRIMMING", "MAPPING", "DEDUP"])
        ):
            subworkflows.remove("QC")

        postprocess = [x for x in wfs if x in post]
        if len(postprocess) == 0 or postprocess[0] == "":
            postprocess = []
    else:
        log.error("NO WORKFLOWS DEFINED, NOTHING TO DO!")
        sys.exit()

    if preprocess:
        try:
            all([config[x] or x == "" for x in preprocess])
        except KeyError:
            log.warning(
                logid
                + "Not all required preprocessing steps have configuration in the config file"
            )

    if subworkflows:
        try:
            all([config[x] or x == "TRIMMING" or x == "" for x in subworkflows])
        except KeyError:
            log.warning(
                logid
                + "Not all required subworkflows have configuration in the config file"
            )

    if postprocess:
        try:
            all([config[x] or x == "" for x in postprocess])
        except KeyError:
            log.warning(
                logid
                + "Not all required postprocessing steps have configuration in the config file"
            )

    log.debug(logid + "WORKFLOWS: " + str([preprocess, subworkflows, postprocess]))

    return [preprocess, subworkflows, postprocess]


@check_run
def nf_make_pre(
    subwork,
    config,
    samples,
    conditions,
    subdir,
    loglevel,
    state="",
    subname=None,
    combinations=None,
):
    logid = scriptname + ".Workflows_nf_make_pre: "
    log.info(
        logid
        + f"PREPROCESSING: {subwork} and SAMPLES: {samples} for COMBINATIONS:{combinations}"
    )

    jobs = list()
    state = "pre_"
    condapath = re.compile(r'conda\s+"')
    includepath = re.compile(r'include:\s+"')
    logfix = re.compile(r'loglevel="INFO"')

    if combinations:
        combname = get_combo_name(combinations)

        for condition in combname:
            worklist = combname[condition]["works"]
            envlist = combname[condition]["envs"]
            add = list()

            nfi = os.path.abspath(os.path.join(workflowpath, "header.nf"))
            with open(nfi, "r") as nf:
                for line in nf.readlines():
                    line = re.sub(logfix, "loglevel='" + loglevel + "'", line)
                    line = re.sub(condapath, 'conda "' + envpath, line)
                    if "include: " in line:
                        line = fixinclude(
                            line, loglevel, condapath, envpath, workflowpath, logfix
                        )
                    add.append(line)
                add.append("\n\n")

            for i in range(len(worklist)):
                log.debug(
                    logid + " SUBLISTS: " + str(worklist[i]) + "\t" + str(envlist[i])
                )
                works = worklist[i].split("-")
                envs = envlist[i].split("-")
                flowlist = list()
                tp = list()
                subjobs = list()
                subconf = NestedDefaultDict()

                for j in range(len(works)):
                    listoftools, listofconfigs = create_subworkflow(
                        config, works[j], [condition], envs
                    )

                    if listoftools is None:
                        log.warning(
                            logid
                            + "No entry fits condition "
                            + str(condition)
                            + " for processing step "
                            + str(works[j])
                        )
                        return None

                    sconf = listofconfigs[0]
                    for a in range(0, len(listoftools)):
                        toolenv, toolbin = map(str, listoftools[a])

                        if toolenv != envs[j] or toolbin is None:
                            continue

                        subsamples = get_samples(sconf)
                        sconf[works[j] + "ENV"] = toolenv
                        sconf[works[j] + "BIN"] = toolbin
                        subconf.merge(sconf)

                        subconf[works[j]] = add_to_innermost_key_by_list(
                            subconf[works[j]],
                            subDict(config[works[j]], condition)[toolenv],
                            condition,
                        )

                        subname = toolenv + ".nf"
                        log.debug(
                            logid
                            + str(works[j])
                            + ": "
                            + str([toolenv, subname, condition, subconf])
                        )

                        if works[j] == "QC":
                            subname = toolenv + "_raw.nf"
                            flowlist.append("QC_RAW")
                            flowlist.append("MULTIQC")

                        nfi = os.path.abspath(os.path.join(workflowpath, subname))
                        with open(nfi, "r") as nf:
                            for line in nf.readlines():
                                line = re.sub(
                                    logfix, "loglevel='" + loglevel + "'", line
                                )
                                line = re.sub(condapath, 'conda "' + envpath, line)
                                if "include: " in line:
                                    line = fixinclude(
                                        line,
                                        loglevel,
                                        condapath,
                                        envpath,
                                        workflowpath,
                                        logfix,
                                    )
                                subjobs.append(line)
                            subjobs.append("\n\n")

                        if works[j] == "QC":
                            nfi = os.path.abspath(
                                os.path.join(workflowpath, "multiqc.nf")
                            )
                            with open(nfi, "r") as nf:
                                for line in nf.readlines():
                                    line = re.sub(
                                        logfix, "loglevel='" + loglevel + "'", line
                                    )
                                    line = re.sub(condapath, 'conda "' + envpath, line)
                                    if "include: " in line:
                                        line = fixinclude(
                                            line,
                                            loglevel,
                                            condapath,
                                            envpath,
                                            workflowpath,
                                            logfix,
                                        )
                                    subjobs.append(line)
                                subjobs.append("\n\n")

                        # workflow merger
                        subjobs.append("\n\n" + "workflow {\n")
                        if "MULTIQC" in flowlist:
                            subjobs.append(" " * 4 + "QC_RAW(dummy)\n")
                            subjobs.append(
                                " " * 4 + "MULTIQC(QC_RAW.out.qc.collect())\n"
                            )
                        subjobs.append("\n}\n")

                        # nfi = os.path.abspath(os.path.join('MONSDA', 'workflows', 'footer.nf'))
                        # with open(nfi, 'r') as nf:
                        #    for line in nf.readlines():
                        #        line = re.sub(logfix, 'loglevel=\''+loglevel+'\'', line)
                        #        line = re.sub(condapath, 'conda: "' + envpath, line)
                        #        subjobs.append(line)
                        #    subjobs.append('\n\n')

                        nfo = os.path.abspath(
                            os.path.join(
                                subdir,
                                "_".join(
                                    [
                                        "_".join(condition),
                                        state + works[j],
                                        toolenv,
                                        "subflow.nf",
                                    ]
                                ),
                            )
                        )
                        if os.path.exists(nfo):
                            os.rename(nfo, nfo + ".bak")
                        with open(nfo, "a") as nfout:
                            nfout.write("".join(subjobs))
                            nfout.write("\n\n")

                        confo = os.path.abspath(
                            os.path.join(
                                subdir,
                                "_".join(
                                    [
                                        "_".join(condition),
                                        state + works[j],
                                        toolenv,
                                        "subconfig.json",
                                    ]
                                ),
                            )
                        )
                        if os.path.exists(confo):
                            os.rename(confo, confo + ".bak")
                        with open(confo, "a") as confout:
                            json.dump(subconf, confout)

                        tp = nf_tool_params(
                            subsamples[0],
                            None,
                            subconf,
                            works[j],
                            toolenv,
                            toolbin,
                            None,
                            condition,
                        )

                        tpl = " ".join(tp)
                        combi = list((str(envlist[i]), ""))
                        para = nf_fetch_params(confo, condition, combi)

                        jobs.append([nfo, confo, tpl, para])

    else:
        for condition in conditions:
            flowlist = list()
            subjobs = list()
            subconf = NestedDefaultDict()
            tp = list()

            log.debug(logid + "PREPARING " + str(subwork) + " " + str(condition))

            nfi = os.path.abspath(os.path.join(workflowpath, "header.nf"))
            with open(nfi, "r") as nf:
                for line in nf.readlines():
                    line = re.sub(logfix, "loglevel='" + loglevel + "'", line)
                    line = re.sub(condapath, 'conda "' + envpath, line)
                    if "include: " in line:
                        line = fixinclude(
                            line, loglevel, condapath, envpath, workflowpath, logfix
                        )
                    subjobs.append(line)
                subjobs.append("\n\n")

            listoftools, listofconfigs = create_subworkflow(
                config, subwork, [condition]
            )
            if listoftools is None:
                log.warning(
                    logid
                    + "No entry fits condition "
                    + str(condition)
                    + " for processing step "
                    + str(subwork)
                )
                return None

            sconf = listofconfigs[0]
            if subwork == "QC":
                subsamples = get_samples(sconf)
            elif subwork == "FETCH":
                subsamples = download_samples(sconf)
            elif subwork == "BASECALL":
                subsamples = basecall_samples(sconf)
            log.debug(logid + f"Running {subwork} for SAMPLES {subsamples}")

            for i in range(0, len(listoftools)):
                toolenv, toolbin = map(str, listoftools[i])
                if toolenv is None or toolbin is None:
                    continue
                sconf[subwork + "ENV"] = toolenv
                sconf[subwork + "BIN"] = toolbin
                subconf.merge(sconf)

                subconf[subwork] = add_to_innermost_key_by_list(
                    subconf[subwork],
                    subDict(config[subwork], condition)[toolenv],
                    condition,
                )

                subname = toolenv + ".nf"
                log.debug(
                    logid
                    + str(subwork)
                    + ": "
                    + str([toolenv, subname, condition, subconf])
                )

                if subwork == "QC":
                    subname = toolenv + "_raw.nf"
                    flowlist.append("QC_RAW")
                elif subwork == "FETCH":
                    subname = toolenv + ".nf"
                    flowlist.append("FETCH")
                elif subwork == "BASECALL":
                    sys.exit("BASECALLING NOT YET IMPLEMENTED FOR NEXTFLOW")
                    subname = toolenv + ".nf"
                    flowlist.append("BASECALL")

                nfi = os.path.abspath(os.path.join(workflowpath, subname))
                with open(nfi, "r") as nf:
                    for line in nf.readlines():
                        line = re.sub(condapath, 'conda "' + envpath, line)
                        if "include: " in line:
                            line = fixinclude(
                                line, loglevel, condapath, envpath, workflowpath, logfix
                            )
                        subjobs.append(line)
                    subjobs.append("\n\n")

                tp.append(
                    nf_tool_params(
                        subsamples[0],
                        None,
                        sconf,
                        subwork,
                        toolenv,
                        toolbin,
                        None,
                        condition,
                    )
                )
            if subwork == "QC":
                flowlist.append("MULTIQC")
                nfi = os.path.abspath(os.path.join(workflowpath, "multiqc.nf"))
                with open(nfi, "r") as nf:
                    for line in nf.readlines():
                        line = re.sub(condapath, 'conda "' + envpath, line)
                        if "include: " in line:
                            line = fixinclude(
                                line, loglevel, condapath, envpath, workflowpath, logfix
                            )
                        subjobs.append(line)
                    subjobs.append("\n\n")

            # workflow merger
            log.debug("FLOWLIST: " + str(flowlist))

            subjobs.append("\n\n" + "workflow {\n")
            for w in ["QC_RAW", "FETCH", "BASECALL"]:
                if w in flowlist:
                    subjobs.append(" " * 4 + w + "(dummy)\n")
            if "MULTIQC" in flowlist:
                subjobs.append(" " * 4 + "MULTIQC(QC_RAW.out.qc.collect())\n")
            subjobs.append("}\n\n")

            nfo = os.path.abspath(
                os.path.join(
                    subdir,
                    "_".join(
                        ["_".join(condition), state + subwork, toolenv, "subflow.nf"]
                    ),
                )
            )
            if os.path.exists(nfo):
                os.rename(nfo, nfo + ".bak")
            with open(nfo, "w") as nfout:
                nfout.write("".join(subjobs))

            confo = os.path.abspath(
                os.path.join(
                    subdir,
                    "_".join(
                        [
                            "_".join(condition),
                            state + subwork,
                            toolenv,
                            "subconfig.json",
                        ]
                    ),
                )
            )
            if os.path.exists(confo):
                os.rename(confo, confo + ".bak")
            with open(confo, "w") as confout:
                json.dump(subconf, confout)

            tpl = " ".join(tp)
            combi = None
            para = nf_fetch_params(confo, condition, combi)

            jobs.append([nfo, confo, tpl, para])

    return jobs


@check_run
def nf_make_sub(
    subworkflows,
    config,
    samples,
    conditions,
    subdir,
    loglevel,
    subname=None,
    combinations=None,
):
    logid = scriptname + ".Workflows_nf_make_sub: "

    log.info(logid + f"STARTING PROCESSING FOR {conditions} and SAMPLES: {samples}")
    jobs = list()
    condapath = re.compile(r'conda\s+"')
    includepath = re.compile(r'include:\s+"')
    logfix = re.compile(r'loglevel="INFO"')

    if combinations:
        combname = get_combo_name(combinations)

        for condition in combname:
            worklist = combname[condition]["works"]
            envlist = combname[condition]["envs"]
            add = list()

            nfi = os.path.abspath(os.path.join(workflowpath, "header.nf"))
            with open(nfi, "r") as nf:
                for line in nf.readlines():
                    line = re.sub(logfix, "loglevel='" + loglevel + "'", line)
                    line = re.sub(condapath, 'conda "' + envpath, line)
                    if "include: " in line:
                        line = fixinclude(
                            line, loglevel, condapath, envpath, workflowpath, logfix
                        )
                    add.append(line)
                add.append("\n\n")

            for i in range(len(worklist)):
                log.debug(
                    logid + " SUBLISTS: " + str(worklist[i]) + "\t" + str(envlist[i])
                )
                works = worklist[i].split("-")
                envs = envlist[i].split("-")
                flowlist = list()
                tp = list()
                subjobs = list()
                subconf = NestedDefaultDict()
                deduptool = None

                for j in range(len(works)):
                    listoftools, listofconfigs = create_subworkflow(
                        config, works[j], [condition], envs
                    )

                    if listoftools is None:
                        log.warning(
                            logid
                            + "No entry fits condition "
                            + str(condition)
                            + " for processing step "
                            + str(works[j])
                        )
                        return None

                    sconf = listofconfigs[0]
                    for a in range(0, len(listoftools)):
                        toolenv, toolbin = map(str, listoftools[a])

                        if toolenv != envs[j] or toolbin is None:
                            continue

                        subsamples = get_samples(sconf)
                        sconf[works[j] + "ENV"] = toolenv
                        sconf[works[j] + "BIN"] = toolbin
                        subconf.merge(sconf)

                        subconf[works[j]] = add_to_innermost_key_by_list(
                            subconf[works[j]],
                            subDict(config[works[j]], condition)[toolenv],
                            condition,
                        )

                        subname = toolenv + ".nf"
                        log.debug(
                            logid
                            + str(works[j])
                            + ": "
                            + str([toolenv, subname, condition, subconf])
                        )

                        if works[j] == "QC":
                            if "TRIMMING" in works:
                                subname = toolenv + "_trim.nf"
                                flowlist.append("QC_RAW")
                                flowlist.append("TRIMMING")
                                flowlist.append("QC_TRIMMING")
                                if "MAPPING" in works:
                                    subname = toolenv + "_trim_map.nf"
                                    flowlist.append("QC_MAP")
                            else:
                                if "DEDUP" in subworkflows:
                                    subname = toolenv + "_dedup.nf"
                                    flowlist.append("QC_RAW")
                                    if toolenv == "umitools":
                                        flowlist.append("DEDUPEXTRACT")
                                    if "MAPPING" in works:
                                        subname = toolenv + "_dedup_map.nf"
                                        flowlist.append("QC_MAP")
                                else:
                                    subname = toolenv + "_raw.nf"
                                    flowlist.append("QC_RAW")
                                    if "MAPPING" in works:
                                        subname = toolenv + ".nf"
                                        flowlist.append("QC_MAP")

                            subname = toolenv + ".nf"
                            flowlist.append("QC_MAPPING")

                        if works[j] == "TRIMMING" and "TRIMMING" not in flowlist:
                            subname = toolenv + ".nf"
                            flowlist.append("TRIMMING")

                        if works[j] == "DEDUP":
                            if toolenv == "umitools":
                                flowlist.append("PREDEDUP")
                                subconf["PREDEDUP"] = "enabled"
                                if "QC" in flowlist:
                                    flowlist.append("QC_DEDUP")
                                subname = toolenv + ".nf"
                                nfi = os.path.abspath(
                                    os.path.join(workflowpath, subname)
                                )
                                with open(nfi, "r") as nf:
                                    for line in nf.readlines():
                                        line = re.sub(
                                            condapath, 'conda "' + envpath, line
                                        )
                                        if "include: " in line:
                                            line = fixinclude(
                                                line,
                                                loglevel,
                                                condapath,
                                                envpath,
                                                workflowpath,
                                                logfix,
                                            )
                                        subjobs.append(line)
                                    subjobs.append("\n\n")

                                tp.append(
                                    nf_tool_params(
                                        subsamples[0],
                                        None,
                                        sconf,
                                        works[j],
                                        toolenv,
                                        toolbin,
                                        None,
                                        condition,
                                    )
                                )
                            subname = toolenv + "_dedup.nf"

                        nfi = os.path.abspath(os.path.join(workflowpath, subname))
                        with open(nfi, "r") as nf:
                            for line in nf.readlines():
                                line = re.sub(condapath, 'conda "' + envpath, line)
                                if "include: " in line:
                                    line = fixinclude(
                                        line,
                                        loglevel,
                                        condapath,
                                        envpath,
                                        workflowpath,
                                        logfix,
                                    )
                                subjobs.append(line)
                            subjobs.append("\n\n")

                        tp.append(
                            nf_tool_params(
                                subsamples[0],
                                None,
                                sconf,
                                works[j],
                                toolenv,
                                toolbin,
                                None,
                                condition,
                            )
                        )

                if "MAPPING" in works:
                    if "QC" not in works:
                        log.debug(logid + "Mapping without QC!")
                    if "TRIMMING" not in works:
                        log.debug(logid + "Simulated read trimming only!")
                        flowlist.append("TRIMMING")
                        nfi = os.path.abspath(
                            os.path.join(
                                installpath,
                                "MONSDA",
                                "workflows",
                                "simulatetrim.nf",
                            )
                        )
                        with open(nfi, "r") as nf:
                            for line in nf.readlines():
                                line = re.sub(condapath, 'conda "' + envpath, line)
                                if "include: " in line:
                                    line = fixinclude(
                                        line,
                                        loglevel,
                                        condapath,
                                        envpath,
                                        workflowpath,
                                        logfix,
                                    )
                                subjobs.append(line)
                            subjobs.append("\n\n")

                    flowlist.append("MAPPING")

                    nfi = os.path.abspath(os.path.join(workflowpath, "mapping.nf"))
                    with open(nfi, "r") as nf:
                        for line in nf.readlines():
                            line = re.sub(condapath, 'conda "' + envpath, line)
                            if "include: " in line:
                                line = fixinclude(
                                    line,
                                    loglevel,
                                    condapath,
                                    envpath,
                                    workflowpath,
                                    logfix,
                                )
                            subjobs.append(line)
                        subjobs.append("\n\n")

                    if "DEDUP" in works:
                        flowlist.append("DEDUPBAM")

                if "QC" in works:
                    flowlist.append("MULTIQC")
                    nfi = os.path.abspath(os.path.join(workflowpath, "multiqc.nf"))
                    with open(nfi, "r") as nf:
                        for line in nf.readlines():
                            line = re.sub(condapath, 'conda "' + envpath, line)
                            if "include: " in line:
                                line = fixinclude(
                                    line,
                                    loglevel,
                                    condapath,
                                    envpath,
                                    workflowpath,
                                    logfix,
                                )
                            subjobs.append(line)
                        subjobs.append("\n\n")

                # workflow merger
                log.debug("FLOWLIST: " + str(flowlist))

                subjobs.append("\n\n" + "workflow {\n")
                for w in [
                    "QC_RAW",
                    "PREDEDUP",
                    "QC_DEDUP",
                    "TRIMMING",
                    "QC_TRIMMING",
                    "MAPPING",
                    "DEDUPBAM",
                    "QC_MAPPING",
                    "MULTIQC",
                ]:
                    if w in flowlist:
                        if w == "QC_RAW":
                            subjobs.append(" " * 4 + w + "(dummy)\n")
                        elif w == "PREDEDUP":
                            subjobs.append(" " * 4 + "DEDUPEXTRACT" + "(dummy)\n")
                        elif w == "QC_DEDUP":
                            subjobs.append(" " * 4 + w + "(DEDUPEXTRACT.out.extract)\n")
                        elif w == "TRIMMING":
                            if "PREDEDUP" in flowlist:
                                subjobs.append(
                                    " " * 4
                                    + "TRIMMING"
                                    + "(DEDUPEXTRACT.out.extract)\n"
                                )
                            else:
                                subjobs.append(" " * 4 + "TRIMMING" + "(dummy)\n")
                        elif w == "QC_TRIMMING":
                            subjobs.append(" " * 4 + w + "(TRIMMING.out.trimmed)\n")
                        elif w == "MAPPING":
                            subjobs.append(" " * 4 + w + "(TRIMMING.out.trimmed)\n")
                            subjobs.append(
                                " " * 4 + "POSTMAPPING(MAPPING.out.mapped)\n"
                            )
                        elif w == "DEDUPBAM":
                            subjobs.append(
                                " " * 4
                                + w
                                + "(POSTMAPPING.out.postmap, POSTMAPPING.out.postbai, POSTMAPPING.out.postmapuni, POSTMAPPING.out.postunibai)\n"
                            )
                        elif w == "QC_MAPPING":
                            if "DEDUPBAM" in flowlist:
                                subjobs.append(
                                    " " * 4
                                    + w
                                    + "(POSTMAPPING.out.postmap.concat(POSTMAPPING.out.postmapuni.concat(DEDUPBAM.out.dedup)))\n"
                                )
                            else:
                                subjobs.append(
                                    " " * 4
                                    + w
                                    + "(POSTMAPPING.out.postmap.concat(POSTMAPPING.out.postmapuni))\n"
                                )
                        elif w == "MULTIQC":
                            if "DEDUPBAM" in flowlist and "QC_TRIMMING" in flowlist:
                                subjobs.append(
                                    " " * 4
                                    + w
                                    + "(QC_RAW.out.qc.concat(QC_TRIMMING.out.qc.concat(QC_MAPPING.out.qc.concat(MAPPING.out.logs))).collect())\n"
                                )
                            elif (
                                "DEDUPBAM" in flowlist and "QC_TRIMMING" not in flowlist
                            ):
                                subjobs.append(
                                    " " * 4
                                    + w
                                    + "(QC_RAW.out.qc.concat(QC_MAPPING.out.qc.concat(MAPPING.out.logs)).collect())\n"
                                )
                            elif "MAPPING" in flowlist and "QC_TRIMMING" in flowlist:
                                subjobs.append(
                                    " " * 4
                                    + w
                                    + "(QC_RAW.out.qc.concat(QC_TRIMMING.out.qc.concat(QC_MAPPING.out.qc.concat(POSTMAPPING.out.postmapuni))).collect())\n"
                                )
                            elif (
                                "MAPPING" in flowlist and "QC_TRIMMING" not in flowlist
                            ):
                                subjobs.append(
                                    " " * 4
                                    + w
                                    + "(QC_RAW.out.qc.concat(QC_MAPPING.out.qc.concat(POSTMAPPING.out.postmapuni)).collect())\n"
                                )
                            elif "TRIMMING" in flowlist and "QC_TRIMMING" in flowlist:
                                subjobs.append(
                                    " " * 4
                                    + w
                                    + "(QC_RAW.out.qc.concat(QC_TRIMMING.out.qc).collect())\n"
                                )
                            else:
                                subjobs.append(
                                    " " * 4 + w + "(QC_RAW.out.qc.collect())\n"
                                )
                        else:
                            subjobs.append(" " * 4 + w + "(dummy)\n")
                subjobs.append("}\n\n")

                # Append footer and write out subflow and subconf per condition
                # nfi = os.path.abspath(os.path.join(installpath, 'MONSDA', 'workflows', 'footer.nf'))
                # with open(nfi,'r') as nf:
                #    for line in nf.readlines():
                #        line = re.sub(condapath,'conda  "../', line)
                #        subjobs.append(line)
                #    subjobs.append('\n\n')

                nfo = os.path.abspath(
                    os.path.join(
                        subdir,
                        "_".join(["_".join(condition), envlist[i], "subflow.nf"]),
                    )
                )
                if os.path.exists(nfo):
                    os.rename(nfo, nfo + ".bak")
                with open(nfo, "w") as nfout:
                    nfout.write("".join(add))
                    nfout.write("".join(subjobs))

                confo = os.path.abspath(
                    os.path.join(
                        subdir,
                        "_".join(["_".join(condition), envlist[i], "subconfig.json"]),
                    )
                )
                if os.path.exists(confo):
                    os.rename(confo, confo + ".bak")
                with open(confo, "w") as confout:
                    json.dump(subconf, confout)

                tpl = " ".join(tp)
                combi = list((str(envlist[i]), ""))
                para = nf_fetch_params(confo, condition, combi)

                jobs.append([nfo, confo, tpl, para])

    else:
        for condition in conditions:
            flowlist = list()
            subjobs = list()
            subconf = NestedDefaultDict()
            tp = list()

            for subwork in subworkflows:
                log.debug(logid + "PREPARING " + str(subwork) + " " + str(condition))
                listoftools, listofconfigs = create_subworkflow(
                    config, subwork, [condition]
                )
                if listoftools is None:
                    log.warning(
                        logid
                        + "No entry fits condition "
                        + str(condition)
                        + " for processing step "
                        + str(subwork)
                    )
                    return None

                sconf = listofconfigs[0]
                subsamples = get_samples(sconf)
                for i in range(0, len(listoftools)):
                    toolenv, toolbin = map(str, listoftools[i])
                    if toolenv is None or toolbin is None:
                        continue
                    sconf[subwork + "ENV"] = toolenv
                    sconf[subwork + "BIN"] = toolbin
                    subconf.merge(sconf)

                    subconf[subwork] = add_to_innermost_key_by_list(
                        subconf[subwork],
                        subDict(config[subwork], condition)[toolenv],
                        condition,
                    )

                    subname = toolenv + ".nf"
                    log.debug(
                        logid
                        + str(subwork)
                        + ": "
                        + str([toolenv, subname, condition, subconf])
                    )

                    if subwork == "QC":
                        if "TRIMMING" in subworkflows:
                            if "DEDUP" in subworkflows:
                                subname = toolenv + "_dedup_trim.nf"
                                flowlist.append("QC_RAW")
                                flowlist.append("DEDUP_TRIM")
                                flowlist.append("QC_DEDUP_TRIM")
                                if "MAPPING" in subworkflows:
                                    subname = toolenv + "_dedup_trim_map.nf"
                                    flowlist.append("QC_MAP")
                            else:
                                subname = toolenv + "_trim.nf"
                                flowlist.append("QC_RAW")
                                flowlist.append("TRIMMING")
                                flowlist.append("QC_TRIMMING")
                                if "MAPPING" in subworkflows:
                                    subname = toolenv + "_trim_map.nf"
                                    flowlist.append("QC_MAP")
                        else:
                            if "DEDUP" in subworkflows:
                                subname = toolenv + "_dedup.nf"
                                flowlist.append("QC_RAW")
                                flowlist.append("DEDUP")
                                flowlist.append("QC_DEDUP")
                                if "MAPPING" in subworkflows:
                                    subname = toolenv + "_dedup_map.nf"
                                    flowlist.append("QC_MAP")
                            else:
                                subname = toolenv + "_raw.nf"
                                flowlist.append("QC_RAW")
                                if "MAPPING" in subworkflows:
                                    subname = toolenv + ".nf"
                                    flowlist.append("QC_MAP")

                        subname = toolenv + ".nf"
                        flowlist.append("QC_MAPPING")

                    # Picard tools can be extended here
                    if subwork == "DEDUP" and toolenv == "picard":
                        subname = toolenv + "_dedup.nf"

                    nfi = os.path.abspath(os.path.join(workflowpath, subname))
                    with open(nfi, "r") as nf:
                        for line in nf.readlines():
                            line = re.sub(condapath, 'conda "' + envpath, line)
                            if "include: " in line:
                                line = fixinclude(
                                    line,
                                    loglevel,
                                    condapath,
                                    envpath,
                                    workflowpath,
                                    logfix,
                                )
                            subjobs.append(line)
                        subjobs.append("\n\n")

                    tp.append(
                        nf_tool_params(
                            samples[0],
                            None,
                            sconf,
                            subwork,
                            toolenv,
                            toolbin,
                            None,
                            condition,
                        )
                    )

            if "MAPPING" in subworkflows:
                if "QC" not in subworkflows:
                    log.debug(logid + "Mapping without QC!")
                if "TRIMMING" not in subworkflows:
                    log.debug(logid + "Simulated read trimming only!")
                    flowlist.append("TRIMMING")
                    nfi = os.path.abspath(os.path.join(workflowpath, "simulatetrim.nf"))
                    with open(nfi, "r") as nf:
                        for line in nf.readlines():
                            line = re.sub(condapath, 'conda "' + envpath, line)
                            if "include: " in line:
                                line = fixinclude(
                                    line,
                                    loglevel,
                                    condapath,
                                    envpath,
                                    workflowpath,
                                    logfix,
                                )
                            subjobs.append(line)
                        subjobs.append("\n\n")

                flowlist.append("MAPPING")

                nfi = os.path.abspath(os.path.join(workflowpath, "mapping.nf"))
                with open(nfi, "r") as nf:
                    for line in nf.readlines():
                        line = re.sub(condapath, 'conda "' + envpath, line)
                        if "include: " in line:
                            line = fixinclude(
                                line, loglevel, condapath, envpath, workflowpath, logfix
                            )
                        subjobs.append(line)
                    subjobs.append("\n\n")

            if "QC" in subworkflows:
                flowlist.append("MULTIQC")
                nfi = os.path.abspath(os.path.join(workflowpath, "multiqc.nf"))
                with open(nfi, "r") as nf:
                    for line in nf.readlines():
                        line = re.sub(condapath, 'conda "' + envpath, line)
                        if "include: " in line:
                            line = fixinclude(
                                line, loglevel, condapath, envpath, workflowpath, logfix
                            )
                        subjobs.append(line)
                    subjobs.append("\n\n")

            # workflow merger
            log.debug("FLOWLIST: " + str(flowlist))

            subjobs.append("\n\n" + "workflow {\n")
            for w in [
                "QC_RAW",
                "TRIMMING",
                "QC_TRIMMING",
                "MAPPING",
                "QC_MAPPING",
                "MULTIQC",
            ]:  # So far DEDUP is missing
                if w in flowlist:
                    if w == "QC_RAW":
                        subjobs.append(" " * 4 + w + "(dummy)\n")
                    elif w == "TRIMMING":
                        if "QC_RAW" not in flowlist:
                            subjobs.append(" " * 4 + "TRIMMING" + "()\n")
                        else:
                            subjobs.append(" " * 4 + "TRIMMING" + "(QC_RAW.out.qc)\n")
                    elif w == "QC_TRIMMING":
                        subjobs.append(" " * 4 + w + "(TRIMMING.out.trimmed)\n")
                    elif w == "MAPPING":
                        subjobs.append(" " * 4 + w + "(TRIMMING.out.trimmed)\n")
                        subjobs.append(" " * 4 + "POSTMAPPING(MAPPING.out.mapped)\n")
                    elif w == "QC_MAPPING":
                        subjobs.append(" " * 4 + w + "(POSTMAPPING.out.postmapuni)\n")
                    elif w == "MULTIQC":
                        if "MAPPING" in flowlist:
                            subjobs.append(" " * 4 + w + "(QC_MAPPING.out.qc)\n")
                        elif "TRIMMING" in flowlist:
                            subjobs.append(" " * 4 + w + "(QC_TRIMMING.out.qc)\n")
                        else:
                            subjobs.append(" " * 4 + w + "(QC_RAW.out.qc)\n")
                    else:
                        subjobs.append(" " * 4 + w + "(dummy)\n")
            subjobs.append("}\n\n")

            # SKIP as nextflow exits too soon otherwise: Append footer and write out subflow and subconf per condition
            # nfi = os.path.abspath(os.path.join(installpath, 'MONSDA', 'workflows', 'footer.nf'))
            # with open(nfi,'r') as nf:
            #    for line in nf.readlines():
            #        line = re.sub(condapath,'conda  "../', line)
            #        subjobs.append(line)
            #    subjobs.append('\n\n')

            nfo = os.path.abspath(
                os.path.join(subdir, "_".join(["_".join(condition), "subflow.nf"]))
            )
            if os.path.exists(nfo):
                os.rename(nfo, nfo + ".bak")
            with open(nfo, "w") as nfout:
                nfout.write("".join(subjobs))

            confo = os.path.abspath(
                os.path.join(subdir, "_".join(["_".join(condition), "subconfig.json"]))
            )
            if os.path.exists(confo):
                os.rename(confo, confo + ".bak")
            with open(confo, "w") as confout:
                json.dump(subconf, confout)

            tpl = " ".join(tp)
            # combi = list((str(envlist[i]), ""))
            combi = None
            para = nf_fetch_params(confo, condition, combi)

            jobs.append([nfo, confo, tpl, para])

    return jobs


@check_run
def nf_make_post(
    postworkflow,
    config,
    samples,
    conditions,
    subdir,
    loglevel,
    subname=None,
    combinations=None,
):
    logid = scriptname + ".Workflows_nf_make_sub: "

    if "PEAKS" in config and "PEAKS" in postworkflow:
        CLIP = checkclip(samples, config)
        log.debug(logid + "Running Peak finding for " + CLIP + " protocol")

    for condition in conditions:
        subconf = NestedDefaultDict()
        for subwork in postworkflow:
            if any(subwork == x for x in ["DE", "DEU", "DAS"]):
                continue
            log.debug(
                logid + "POSTPROCESS: " + str(subwork) + " CONDITION: " + str(condition)
            )
            listoftools, listofconfigs = create_subworkflow(
                config, subwork, [condition], stage="POST"
            )
            log.debug(logid + str([listoftools, listofconfigs]))
            if listoftools is None:
                log.warning(
                    logid
                    + "No entry fits condition "
                    + str(condition)
                    + " for postprocessing step "
                    + str(subwork)
                )
                continue

            for i in range(0, len(listoftools)):
                toolenv, toolbin = map(str, listoftools[i])
                subconf.update(listofconfigs[i])
                subname = toolenv + ".nf"
                subsamples = list(set(sampleslong(subconf)))
                log.debug(
                    logid
                    + "POSTPROCESS: "
                    + str([toolenv, subname, condition, subsamples, subconf])
                )

                nf = os.path.abspath(os.path.join(workflowpath, "header.nf"))
                nfo = os.path.abspath(
                    os.path.join(
                        subdir,
                        "_".join(["_".join(condition), subwork, toolbin, "subflow.nf"]),
                    )
                )
                if os.path.exists(nfo):
                    os.rename(nfo, nfo + ".bak")
                with open(nfo, "a") as nfout:
                    with open(nf, "r") as nfi:
                        for line in nfi.readlines():
                            nfout.write(line)
                    nfout.write("\n\n")
                nfo = os.path.abspath(os.path.join(workflowpath, subname))
                with open(nfo, "a") as nfout:
                    with open(nf, "r") as nfi:
                        smkout.write(smk.read())
                    smkout.write("\n\n")

                # workflow merger
                with open(smko, "a") as smkout:
                    smkout.write("\n\n" + "workflow {\n    main:\n")
                    for w in [
                        "QC_RAW",
                        "TRIMMING",
                        "QC_TRIMMING",
                        "MAPPING",
                        "QC_MAPPING",
                        "MULTIQC",
                    ]:
                        if w in flowlist:
                            smkout.write("\n    ".w)
                    smkout.write("\n}\n")

                # smkf = os.path.abspath(os.path.join(installpath, 'MONSDA','workflows','footer.nf'))
                # with open(smko, 'a') as smkout:
                #    with open(smkf,'r') as smk:
                #        smkout.write(smk.read())

                confo = os.path.abspath(
                    os.path.join(
                        subdir,
                        "_".join(
                            ["_".join(condition), subwork, toolbin, "subconfig.json"]
                        ),
                    )
                )
                if os.path.exists(confo):
                    os.rename(confo, confo + ".bak")

                with open(confo, "a") as confout:
                    json.dump(subconf, confout)

                params = nf_fetch_params(
                    os.path.abspath(
                        os.path.join(
                            subdir, "_".join(["_".join(condition), "subconfig.json"])
                        ),
                        condition,
                        combination,
                    )
                )
                toolparams = nf_tool_params(
                    subsamples[0], None, subconf, subwork, condition
                )

                jobtorun = "nextflow -log /dev/stderr run {s} -w {d} {rest} {p} {j} {c}".format(
                    t=threads,
                    s=os.path.abspath(
                        os.path.join(
                            subdir,
                            "_".join(
                                ["_".join(condition), subwork, toolbin, "subflow.nf"]
                            ),
                        )
                    ),
                    d=workdir,
                    rest=" ".join(argslist),
                    p=" ".join(
                        "--{!s} {!s}".format(key, val) for (key, val) in params.items()
                    ),
                    j=toolparams,
                    c="--CONDITION " + str.join(os.sep, condition),
                )

                log.info(logid + "RUNNING " + str(jobtorun))
                job = runjob(jobtorun)
                log.debug(logid + "JOB CODE " + str(job))

    # THIS SECTION IS FOR DE, DEU, DAS ANALYSIS, WE USE THE CONDITIONS TO MAKE PAIRWISE COMPARISONS
    for analysis in ["DE", "DEU", "DAS"]:
        if analysis in config and analysis in postworkflow:
            log.info(logid + "STARTING " + analysis + " Analysis...")
            subwork = analysis
            subconf = NestedDefaultDict()
            log.debug(
                logid + "SUBWORK: " + str(subwork) + " CONDITION: " + str(conditions)
            )
            listoftools, listofconfigs = create_subworkflow(
                config, subwork, conditions, stage="POST"
            )

            if listoftools is None:  # or listoftoolscount is None:
                log.error(
                    logid
                    + "No entry fits condition "
                    + str(conditions)
                    + " for postprocessing step "
                    + str(subwork)
                )

            for key in config[subwork]["TOOLS"]:
                log.info(logid + "... with Tool: " + key)
                toolenv = key
                toolbin = config[subwork]["TOOLS"][key]
                # countenv, countbin = map(str, listoftoolscount[0]) #Counting per analysis rule now
                subconf = NestedDefaultDict()
                for i in listofconfigs:
                    i[subwork + "ENV"] = toolenv
                    i[subwork + "BIN"] = toolbin
                    # i['COUNTBIN'] = 'featureCounts'#This is hard coded where needed for now
                    # i['COUNTENV'] = 'countreads'#This is hard coded where needed for now
                for i in range(len(listoftools)):
                    subconf = merge_dicts(subconf, listofconfigs[i])

                # for x in range(0, len(listofconfigscount)): ### muss hier auch noch gefiltert werden?
                #    subconf = merge_dicts(subconf, listofconfigscount[x])
                subname = (
                    toolenv + ".nf"
                    if toolenv != "edger"
                    else toolenv + "_" + subwork + ".nf"
                )
                subsamples = sampleslong(subconf)
                log.debug(
                    logid
                    + "POSTPROCESS: "
                    + str([toolenv, subname, subsamples, subconf])
                )

                smkf = os.path.abspath(os.path.join(workflowpath, "header.nf"))
                smko = os.path.abspath(
                    os.path.join(subdir, "_".join([subwork, toolenv, "subflow.nf"]))
                )
                if os.path.exists(smko):
                    os.rename(smko, smko + ".bak")
                with open(smko, "a") as smkout:
                    with open(smkf, "r") as smk:
                        for line in smk.readlines():
                            line = re.sub(condapath, 'conda  "../', line)
                            if "include: " in line:
                                line = fixinclude(
                                    line,
                                    loglevel,
                                    condapath,
                                    envpath,
                                    workflowpath,
                                    logfix,
                                )
                            smkout.write(line)
                    smkout.write("\n\n")
                smkf = os.path.abspath(os.path.join(workflowpath, subname))
                with open(
                    os.path.abspath(
                        os.path.join(subdir, "_".join([subwork, toolenv, "subflow.nf"]))
                    ),
                    "a",
                ) as smkout:
                    with open(smkf, "r") as smk:
                        for line in smk.readlines():
                            line = re.sub(logfix, "loglevel='" + loglevel + "'", line)
                            line = re.sub(condapath, 'conda:  "' + envpath, line)
                            if "include: " in line:
                                line = fixinclude(
                                    line,
                                    loglevel,
                                    condapath,
                                    envpath,
                                    workflowpath,
                                    logfix,
                                )
                            smkout.write(line)
                    smkout.write("\n")

                # smkf = os.path.abspath(os.path.join(installpath, 'MONSDA','workflows','footer.nf'))
                # with open(smko, 'a') as smkout:
                #    with open(smkf,'r') as smk:
                #        smkout.write(smk.read())

                confo = os.path.abspath(
                    os.path.join(subdir, "_".join([subwork, toolenv, "subconfig.json"]))
                )
                if os.path.exists(confo):
                    os.rename(confo, confo + ".bak")
                with open(confo, "a") as confout:
                    json.dump(subconf, confout)

                params = nf_fetch_params(
                    os.path.abspath(
                        os.path.join(
                            subdir, "_".join(["_".join(condition), "subconfig.json"])
                        )
                    ),
                    condition,
                    combination,
                )
                toolparams = nf_tool_params(
                    subsamples[0], None, subconf, subwork, condition
                )

                jobtorun = "nextflow -log /dev/stderr run {s} -w {d} {rest} {p} {j} {c}".format(
                    t=threads,
                    s=os.path.abspath(
                        os.path.join(
                            subdir,
                            "_".join(
                                ["_".join(condition), subwork, toolbin, "subflow.nf"]
                            ),
                        )
                    ),
                    d=workdir,
                    rest=" ".join(argslist),
                    p=" ".join(
                        "--{!s} {!s}".format(key, val) for (key, val) in params.items()
                    ),
                    j=toolparams,
                    c="--CONDITION " + str.join(os.sep, condition),
                )

                log.info(logid + "RUNNING " + str(jobtorun))
                job = runjob(jobtorun)
                log.debug(logid + "JOB CODE " + str(job))


#
# Workflows.py ends here
