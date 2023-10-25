#!/usr/bin/env python3

import argparse
import os
import shlex
import shutil
import subprocess
import sys
import traceback as tb

from snakemake import load_configfile
from snakemake.utils import min_version

scriptname = os.path.basename(__file__).replace("Run", "").replace(".py", "")

from MONSDA.Logger import makelogdir, setup_logger

from . import _version

__version__ = _version.get_versions()["version"]

# Logging
import datetime

makelogdir("LOGS")
if not os.path.isfile(os.path.abspath("LOGS" + os.sep + scriptname + ".log")):
    open("LOGS/" + scriptname + ".log", "w").close()
else:
    ts = str(
        datetime.datetime.fromtimestamp(
            os.path.getmtime(os.path.abspath("LOGS" + os.sep + scriptname + ".log"))
        ).strftime("%Y%m%d_%H_%M_%S")
    )
    shutil.copy2(
        "LOGS" + os.sep + scriptname + ".log",
        "LOGS" + os.sep + scriptname + "_" + ts + ".log",
    )

log = setup_logger(
    name=scriptname.lower(),
    log_file="LOGS" + os.sep + scriptname + ".log",
    logformat="%(asctime)s %(levelname)-8s %(name)-12s %(message)s",
    datefmt="%m-%d %H:%M",
)
log = setup_logger(
    name=scriptname.lower(),
    log_file="stderr",
    logformat="%(asctime)s %(levelname)-8s %(message)s",
    datefmt="%m-%d %H:%M",
)

import MONSDA.Params as mp
import MONSDA.Utils as mu
import MONSDA.Workflows as mw


# CODE
def parseargs():
    parser = argparse.ArgumentParser(
        description="Modular Organizer of Nextflow and Snakemake driven hts Data Analysis"
    )
    parser.add_argument(
        "-c", "--configfile", type=str, help="Configuration json to read"
    )
    parser.add_argument(
        "-d", "--directory", type=str, default="", help="Working Directory"
    )
    parser.add_argument(
        "-u",
        "--use-conda",
        action="store_true",
        default=True,
        help="Should conda be used, default True",
    )
    parser.add_argument(
        "-l",
        "--unlock",
        action="store_true",
        help="If Snakemake directory is locked you can unlock before processing",
    )
    parser.add_argument(
        "-j",
        "--procs",
        type=int,
        default=1,
        help="Number of parallel processes to start MONSDA with, capped by MAXTHREADS in config!",
    )
    parser.add_argument(
        "--save",
        action="store_true",
        help="Do not actually run jobs, create corresponding text file containing CLI-calls and arguments for manual running instead",
    )
    parser.add_argument(
        "-s",
        "--skeleton",
        action="store_true",
        help="Just create the minimal directory hierarchy as needed",
    )
    parser.add_argument(
        "--snakemake",
        action="store_true",
        default=True,
        help="Wrap around Snakemake, default",
    )
    parser.add_argument(
        "--nextflow", action="store_true", default=False, help="Wrap around Nextflow"
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Cleanup workdir (Nextflow), append -n to see list of files to clean or -f to actually remove those files",
    )
    parser.add_argument(
        "--loglevel",
        type=str,
        default="INFO",
        choices=["WARNING", "ERROR", "INFO", "DEBUG"],
        help="Set log level",
    )
    parser.add_argument(
        "--version",
        action="store_true",
        help="Print version and exit",
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_known_args()


def run_snakemake(
    configfile,
    workdir,
    useconda,
    procs,
    skeleton,
    loglevel,
    save=None,
    unlock=None,
    optionalargs=None,
):
    try:
        logid = scriptname + ".run_snakemake: "
        config = load_configfile(configfile)
        subdir = "SubSnakes"
        mp.create_skeleton(subdir, skeleton)

        argslist = list()
        if useconda:
            argslist.append("--use-conda")
        else:
            log.warning(
                logid
                + "You are not making use of conda, be aware that this will most likely not work for the workflows provided here! To change append the --use-conda option to your commandline call. You can also speed up conda with the --conda-frontend mamba argument and preinstall all conda environments appending the --use-conda and the --create-envs-only arguments and share conda environment locations across runs with the --conda-prefix argument."
            )
        if optionalargs and len(optionalargs) > 0:
            log.debug(logid + "OPTIONALARGS: " + str(optionalargs))
            argslist.extend(optionalargs)
            if "--profile" in optionalargs and "MONSDA/slurm" in optionalargs:
                mu.makeoutdir("LOGS/SLURM")

        threads = (
            min(int(config["MAXTHREADS"]), procs) if "MAXTHREADS" in config else procs
        )
        config["MAXTHREADS"] = threads

        if unlock:
            try:
                pythonversion = (
                    f"python{str(sys.version_info.major)}.{str(sys.version_info.minor)}"
                )
                installpath = os.path.dirname(__file__).replace(
                    os.sep.join(["lib", pythonversion, "site-packages", "MONSDA"]),
                    "share",
                )
            except:
                installpath = os.getcwd()

            workflowpath = os.path.join(installpath, "MONSDA", "workflows")
            log.info(logid + "Unlocking directory")
            snk = os.path.abspath(os.path.join(workflowpath, "unlock.smk"))
            jobtorun = (
                f"snakemake --unlock -j {threads} -s {snk} --configfile {configfile}"
            )

            log.info(logid + "UNLOCKING " + str(jobtorun))
            job = runjob(jobtorun)
            log.debug(logid + "JOB CODE " + str(job))

        # Get processes to work on
        preprocess, subworkflows, postprocess = mw.get_processes(config)
        conditions = mp.get_conditions(config)

        """
        START TO PREPROCESS
        IF WE NEED TO DOWNLOAD FILES WE DO THIS NOW
        """
        if preprocess:
            for proc in [
                x for x in preprocess if config.get(x) and x in ["FETCH", "BASECALL"]
            ]:
                log.info(logid + "Preprocess " + str(proc))
                if proc not in config:
                    log.error(
                        logid
                        + "No configuration with key "
                        + proc
                        + " for file download found. Nothing to do!"
                    )
                mu.makeoutdir("FASTQ")
                mu.makeoutdir("TMP")

                if proc == "FETCH":
                    SAMPLES = mp.download_samples(config)
                    preprocess.remove(proc)
                elif proc == "BASECALL":
                    SAMPLES = mp.basecall_samples(config)
                    preprocess.remove(proc)
                else:
                    continue  # We only want download/basecall here

                log.debug(logid + "PRESAMPLES: " + str(SAMPLES))
                combinations = (
                    mw.get_combo(preprocess, config, conditions)
                    if subworkflows
                    else None
                )

                subwork = proc

                jobs = mw.make_pre(
                    subwork,
                    config,
                    SAMPLES,
                    conditions,
                    subdir,
                    loglevel,
                    "Pre",
                    combinations=combinations,
                )

                jobstorun = list()
                for job in jobs:
                    smko, confo = job
                    rest = " ".join(argslist)
                    jobstorun.append(
                        f"snakemake -j {threads} -s {smko} --configfile {confo} --directory {workdir} --printshellcmds --show-failed-logs {rest}"
                    )

                for job in jobstorun:
                    with open("JOBS" + os.sep + scriptname + ".commands", "a") as j:
                        j.write(job + os.linesep)
                    if not save:
                        log.info(logid + "RUNNING " + str(job))
                        jid = runjob(job)
                        log.debug(logid + "JOB CODE " + str(jid))

        """
        ONCE FILES ARE DOWNLOADED WE CAN START OTHER PREPROCESSING STEPS
        """
        mu.makeoutdir("TMP")
        SAMPLES = mp.get_samples(config)
        log.info(logid + "SAMPLES: " + str(SAMPLES))

        if preprocess:
            log.info(logid + "STARTING PREPROCESSING")
            if "QC" in preprocess and "QC" in config:
                mu.makeoutdir("QC")

                for subwork in preprocess:
                    combinations = (
                        mw.get_combo(subworkflows, config, conditions)
                        if subworkflows
                        else None
                    )
                    jobs = mw.make_pre(
                        subwork,
                        config,
                        SAMPLES,
                        conditions,
                        subdir,
                        loglevel,
                        "Pre",
                        combinations=combinations,
                    )

                    jobstorun = list()

                    for job in jobs:
                        smko, confo = job
                        rest = " ".join(argslist)

                        jobstorun.append(
                            f"snakemake -j {threads} --use-conda -s {smko} --configfile {confo} --directory {workdir} --printshellcmds --show-failed-logs {rest}"
                        )

                    for job in jobstorun:
                        with open("JOBS" + os.sep + scriptname + ".commands", "a") as j:
                            j.write(job + os.linesep)
                        if not save:
                            log.info(logid + "RUNNING " + str(job))
                            jid = runjob(job)
                            log.debug(logid + "JOB CODE " + str(jid))

        else:
            log.warning(
                logid + "No preprocessing workflows defined! Continuing with workflows!"
            )

        """
        END OF PREPROCESSING, START OF PROCESSING
        """
        mu.makeoutdir("TMP")
        if subworkflows:
            combinations = (
                mw.get_combo(subworkflows, config, conditions) if subworkflows else None
            )
            jobs = mw.make_sub(
                subworkflows,
                config,
                SAMPLES,
                conditions,
                subdir,
                loglevel,
                combinations=combinations,
            )

            jobstorun = list()

            for job in jobs:
                smko, confo = job
                rest = " ".join(argslist)

                jobstorun.append(
                    f"snakemake -j {threads} --use-conda -s {smko} --configfile {confo} --directory {workdir} --printshellcmds --show-failed-logs {rest}"
                )

            for job in jobstorun:
                with open("JOBS" + os.sep + scriptname + ".commands", "a") as j:
                    j.write(job + os.linesep)
                    if not save:
                        log.info(logid + "RUNNING " + str(job))
                        jid = runjob(job)
                        log.debug(logid + "JOB CODE " + str(jid))

        else:
            log.warning(
                logid
                + "No Workflows defined! Nothing to do, continuing with postprocessing!"
            )

        """
        END OF PROCESSING, START OF POSTPROCESSING
        """
        mu.makeoutdir("TMP")
        if postprocess:
            for subwork in postprocess:
                SAMPLES = mp.get_samples_postprocess(config, subwork)
                log.info(logid + f"POSTPROCESSING {subwork} with SAMPLES: {SAMPLES}")
                combinations = (
                    mw.get_combo(subworkflows, config, conditions)
                    if subworkflows
                    else None
                )

                jobs = mw.make_post(
                    subwork,
                    config,
                    SAMPLES,
                    conditions,
                    subdir,
                    loglevel,
                    combinations=combinations,
                )
                jobstorun = list()

                for job in jobs:
                    smko, confo = job
                    rest = " ".join(argslist)
                    jobstorun.append(
                        f"snakemake -j {threads} --use-conda -s {smko} --configfile {confo} --directory {workdir} --printshellcmds --show-failed-logs {rest}"
                    )

                for job in jobstorun:
                    with open("JOBS" + os.sep + scriptname + ".commands", "a") as j:
                        j.write(job + os.linesep)
                        if not save:
                            log.info(logid + "RUNNING " + str(job))
                            jid = runjob(job)
                            log.debug(logid + "JOB CODE " + str(jid))

            if (
                any([x in postprocess for x in ["DE", "DEU", "DAS", "DTU"]])
                and not save
            ):
                # SUMMARY RUN
                combinations = (
                    mw.get_combo(subworkflows, config, conditions)
                    if subworkflows
                    else None
                )
                jobs = mw.make_summary(
                    config, subdir, loglevel, combinations=combinations
                )
                jobstorun = list()

                for job in jobs:
                    smko, confo = job
                    rest = " ".join(argslist)
                    jobstorun.append(
                        f"snakemake -j {threads} --use-conda -s {smko} --configfile {confo} --directory {workdir} --printshellcmds --show-failed-logs {rest}"
                    )

                for job in jobstorun:
                    with open("JOBS" + os.sep + scriptname + ".commands", "a") as j:
                        j.write(job + os.linesep)
                        if not save:
                            log.info(logid + "RUNNING " + str(job))
                            jid = runjob(job)
                            log.debug(logid + "JOB CODE " + str(jid))

        else:
            log.warning(logid + "No postprocessing steps defined! Nothing to do!")

        if save:
            log.info(
                f"{logid} All CLI calls have been saved to monsda.commands in directory JOBS"
            )

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error("".join(tbe.format()))


def run_nextflow(
    configfile,
    workdir,
    useconda,
    procs,
    skeleton,
    loglevel,
    save=None,
    clean=None,
    optionalargs=None,
):
    try:
        logid = scriptname + ".run_nextflow: "
        argslist = list()
        if useconda:
            argslist.append("-with-conda")
        else:
            log.warning(
                logid
                + "You are not making use of conda, be aware that this will most likely not work for the workflows provided here! To change append the --use-conda option to your commandline call."
            )
        if optionalargs and len(optionalargs) > 0:
            log.debug(logid + "OPTIONALARGS: " + str(optionalargs))
            argslist.extend(optionalargs)

        if clean:
            log.info(logid + "Cleaning working directory")
            a = " ".join(argslist)
            jobtorun = f"nextflow clean {a}"
            log.info(logid + "CLEANUP " + str(jobtorun))
            job = runjob(jobtorun)
            log.debug(logid + "JOB CODE " + str(job))
            sys.exit()

        config = load_configfile(configfile)
        workdir = os.path.abspath(str.join(os.sep, [workdir, "NextFlowWork"]))
        subdir = "SubFlows"
        mp.create_skeleton(subdir, skeleton)

        if not os.path.exists(os.path.abspath(subdir + os.sep + "bin")):
            os.symlink(
                os.path.abspath(config.get("BINS", "")),
                os.path.abspath(subdir + os.sep + "bin"),
            )

        threads = (
            min(int(config["MAXTHREADS"]), procs) if "MAXTHREADS" in config else procs
        )
        config["MAXTHREADS"] = threads

        # Get processes to work on
        preprocess, subworkflows, postprocess = mw.nf_get_processes(config)
        conditions = mp.get_conditions(config)

        """
        START TO PROCESS
        IF WE NEED TO DOWNLOAD FILES WE DO THIS NOW
        """

        if preprocess:
            for proc in [
                x for x in preprocess if config.get(x) and x in ["FETCH", "BASECALL"]
            ]:
                log.info(logid + "Preprocess " + str(proc))
                if proc not in config:
                    log.error(
                        logid
                        + "No configuration with key "
                        + proc
                        + " for file download found. Nothing to do!"
                    )
                mu.makeoutdir("FASTQ")
                mu.makeoutdir("TMP")

                if proc == "FETCH":
                    if mp.check_samples(config) is True:
                        preprocess.remove(proc)
                        SAMPLES = None
                    else:
                        SAMPLES = mp.download_samples(config)
                        preprocess.remove(proc)
                elif proc == "BASECALL":
                    SAMPLES = mp.basecall_samples(config)
                    preprocess.remove(proc)
                else:
                    continue
                if SAMPLES is None:
                    log.info(
                        logid + f"All SAMPLES already available, skipping FETCH process"
                    )
                    continue

                log.debug(logid + "PRESAMPLES: " + str(SAMPLES))
                combinations = (
                    mw.get_combo(preprocess, config, conditions)
                    if subworkflows
                    else None
                )

                subwork = proc

                jobs = mw.nf_make_pre(
                    subwork,
                    config,
                    SAMPLES,
                    conditions,
                    subdir,
                    loglevel,
                    "Pre",
                    combinations=combinations,
                )

                jobstorun = list()
                for job in jobs:
                    nfo, confo, tp, params = job
                    pars = " ".join(
                        "--{!s} {!s}".format(key, val) for (key, val) in params.items()
                    )
                    rest = " ".join(argslist)

                    jobstorun.append(
                        f"nextflow -log {workdir}/run.log run {nfo} -w {workdir} {rest} {pars} {tp}"
                    )

                for job in jobstorun:
                    with open("JOBS" + os.sep + scriptname + ".commands", "a") as j:
                        j.write(job + os.linesep)
                    if not save:
                        log.info(logid + "RUNNING " + str(job))
                        jid = runjob(job)
                        log.debug(logid + "JOB CODE " + str(jid))

        """
        ONCE FILES ARE DOWNLOAD WE CAN START OTHER PREPROCESSING STEPS
        """
        mu.makeoutdir("TMP")
        SAMPLES = mp.get_samples(config)
        log.info(logid + "SAMPLES: " + str(SAMPLES))
        conditions = mp.get_conditions(config)
        log.info(logid + "CONDITIONS: " + str(conditions))

        if preprocess:
            log.info(logid + "STARTING PREPROCESSING")
            if "QC" in preprocess and "QC" in config:
                mu.makeoutdir("QC")

            for subwork in preprocess:
                combinations = (
                    mw.get_combo(subworkflows, config, conditions)
                    if subworkflows
                    else None
                )
                jobs = mw.nf_make_pre(
                    subwork,
                    config,
                    SAMPLES,
                    conditions,
                    subdir,
                    loglevel,
                    combinations=combinations,
                )

                jobstorun = list()

                for job in jobs:
                    nfo, confo, tp, params = job
                    pars = " ".join(
                        "--{!s} {!s}".format(key, val) for (key, val) in params.items()
                    )
                    rest = " ".join(argslist)

                    jobstorun.append(
                        f"nextflow -log {workdir}/run.log run {nfo} -w {workdir} {rest} {pars} {tp}"
                    )
                    log.debug(f"{logid} JOBS: {jobstorun}")

                for job in jobstorun:
                    with open("JOBS" + os.sep + scriptname + ".commands", "a") as j:
                        j.write(job + os.linesep)
                    if not save:
                        log.info(logid + "RUNNING " + str(job))
                        jid = runjob(job)
                        log.debug(logid + "JOB CODE " + str(jid))

        else:
            log.warning(
                logid + "No preprocessing workflows defined! Continuing with workflows!"
            )

        """
        END OF PREPROCESSING, START OF PROCESSING
        """
        mu.makeoutdir("TMP")
        if subworkflows:
            combinations = (
                mw.get_combo(subworkflows, config, conditions) if subworkflows else None
            )
            jobs = mw.nf_make_sub(
                subworkflows,
                config,
                SAMPLES,
                conditions,
                subdir,
                loglevel,
                combinations=combinations,
            )

            jobstorun = list()

            for job in jobs:
                nfo, confo, tp, params = job
                pars = " ".join(
                    "--{!s} {!s}".format(key, val) for (key, val) in params.items()
                )
                rest = " ".join(argslist)

                jobstorun.append(
                    f"nextflow -log {workdir}/run.log run {nfo} -w {workdir} {rest} {pars} {tp}"
                )

            log.debug(f"{logid} JOBS: {jobstorun}")
            for job in jobstorun:
                with open("JOBS" + os.sep + scriptname + ".commands", "a") as j:
                    j.write(job + os.linesep)
                    if not save:
                        log.info(logid + "RUNNING " + str(job))
                        jid = runjob(job)
                        log.debug(logid + "JOB CODE " + str(jid))

        else:
            log.warning(
                logid
                + "No Workflows defined! Nothing to do, continuing with postprocessing!"
            )

        """
        END OF PROCESSING, START OF POSTPROCESSING
        """
        mu.makeoutdir("TMP")
        if postprocess:
            for subwork in postprocess:
                SAMPLES = mp.get_samples_postprocess(config, subwork)
                log.info(logid + "POSTPROCESSING SAMPLES: " + str(SAMPLES))
                combinations = (
                    mw.get_combo(subworkflows, config, conditions)
                    if subworkflows
                    else None
                )
                log.debug(logid + "POSTPROCESSING WITH COMBOS: " + str(combinations))

                jobs = mw.nf_make_post(
                    subwork,
                    config,
                    SAMPLES,
                    conditions,
                    subdir,
                    loglevel,
                    combinations=combinations,
                )
                jobstorun = list()

                for job in jobs:
                    nfo, confo, tp, params = job
                    pars = " ".join(
                        "--{!s} {!s}".format(key, val) for (key, val) in params.items()
                    )
                    rest = " ".join(argslist)

                    jobstorun.append(
                        f"nextflow -log {workdir}/run.log run {nfo} -w {workdir} {rest} {pars} {tp}"
                    )

                for job in jobstorun:
                    with open("JOBS" + os.sep + scriptname + ".commands", "a") as j:
                        j.write(job + os.linesep)
                        if not save:
                            log.info(logid + "RUNNING " + str(job))
                            jid = runjob(job)
                            log.debug(logid + "JOB CODE " + str(jid))

            if (
                any([x in postprocess for x in ["DE", "DEU", "DAS", "DTU"]])
                and not save
            ):
                # SUMMARY RUN
                combinations = (
                    mw.get_combo(subworkflows, config, conditions)
                    if subworkflows
                    else None
                )
                jobs = mw.nf_make_summary(
                    config, subdir, loglevel, combinations=combinations
                )
                jobstorun = list()

                for job in jobs:
                    nfo, confo, tp, params = job
                    pars = " ".join(
                        "--{!s} {!s}".format(key, val) for (key, val) in params.items()
                    )
                    rest = " ".join(argslist)

                    jobstorun.append(
                        f"nextflow -log {workdir}/run.log run {nfo} -w {workdir} {rest} {pars} {tp}"
                    )

                for job in jobstorun:
                    with open("JOBS" + os.sep + scriptname + ".commands", "a") as j:
                        j.write(job + os.linesep)
                        if not save:
                            log.info(logid + "RUNNING " + str(job))
                            jid = runjob(job)
                            log.debug(logid + "JOB CODE " + str(jid))

        else:
            log.warning(logid + "No postprocessing steps defined! Nothing to do!")

        if save:
            log.info(
                f"{logid} All CLI calls have been saved to monsda.commands in directory JOBS"
            )

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error("".join(tbe.format()))


def runjob(jobtorun):
    try:
        logid = scriptname + ".runjob: "
        # return subprocess.run(jobtorun, shell=True, universal_newlines=True, capture_output=True)  # python >= 3.7
        job = subprocess.Popen(
            jobtorun,
            shell=True,
            universal_newlines=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            bufsize=1,
            close_fds=True,
        )

        while True:
            output, outerr = job.communicate()
            output = str.join("", output).rstrip()
            outerr = str.join("", outerr).rstrip()
            if output == "" and outerr == "" and job.poll() is not None:
                break
            if output and output != "":
                if (
                    any(x in output for x in ["ERROR", "Error", "error", "Exception"])
                    and not "Workflow finished" in output
                    or "Workflow finished" in outerr
                ):
                    if outerr:
                        log.error(
                            logid + "STOPPING: " + str(output) + "\n" + str(outerr)
                        )
                    else:
                        log.error(logid + "STOPPING: " + str(output))
                    log.info("PLEASE CHECK LOG AT LOGS/MONSDA.log")
                    job.kill()
                    sys.exit()
                else:
                    log.info(logid + str(output))
            if outerr and outerr != "":
                if (
                    not "Workflow finished" in outerr
                    and not "Nothing to be done" in outerr
                    and not "Workflow finished" in output
                    and any(
                        x in outerr for x in ["ERROR", "Error", "error", "Exception"]
                    )
                ):
                    log.error(logid + "STOPPING: " + str(outerr))
                    log.info("PLEASE CHECK LOG AT LOGS/MONSDA.log")
                    job.kill()
                    sys.exit()
                else:
                    log.info(logid + str(outerr))
            if job.poll() is not None:
                break

        if job.returncode == 0:
            output, outerr = job.communicate()
            output = str.join("", output).rstrip()
            if output and output != "":
                log.info(logid + "JOB FINISHED: " + output)
            return job.poll()
        else:
            output, outerr = job.communicate()
            output = str.join("", output).rstrip()
            outerr = str.join("", outerr).rstrip()
            if outerr and outerr != "" or output and output != "":
                log.error(logid + "ERROR: " + outerr + output)
                log.info("PLEASE CHECK LOG AT LOGS/MONSDA.log")
            job.kill()
            sys.exit("ERROR SIGNAL: " + str(job.returncode))

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error("".join(tbe.format()))
        sys.exit()


def main():
    logid = scriptname + ".main: "
    try:
        args = parseargs()
        knownargs = args[0]
        optionalargs = args[1:]

        if knownargs.version:
            sys.exit("MONSDA version " + __version__)

        log.setLevel(knownargs.loglevel)

        log.info(
            logid
            + "MONSDA-CLI: "
            + sys.argv[0]
            + " "
            + "{}".format(" ".join([shlex.quote(s) for s in sys.argv[1:]]))
        )

        required_version = load_configfile(knownargs.configfile).get("VERSION")
        if not required_version:
            sys.exit(
                "Can not check version needed, please add VERSION key to config file"
            )

        if mw.monsda_check_version(__version__, required_version) is not True:
            log.error(
                "Version required in config file "
                + str(required_version)
                + " does not match installed version "
                + str(__version__)
                + " Please install correct version before continuing!"
            )
            sys.exit(
                "Version required in config file "
                + str(required_version)
                + " does not match installed version "
                + str(__version__)
                + " Please install correct version before continuing!"
            )
        else:
            log.info("Running MONSDA version " + __version__ + " as configured")

        MIN_PYTHON = (3, 11)
        if sys.version_info < MIN_PYTHON:
            log.error("This script requires Python version >= 3.7")
            sys.exit("This script requires Python version >= 3.7")
        log.info(
            logid + "Running " + scriptname + " on " + str(knownargs.procs) + " cores"
        )
        log.debug(logid + str(log.handlers))

        commandfile = "JOBS" + os.sep + scriptname + ".commands"
        if os.path.exists(commandfile):
            ts = "_" + str(
                datetime.datetime.fromtimestamp(
                    os.path.getmtime(os.path.abspath(commandfile))
                ).strftime("%Y%m%d_%H_%M_%S")
            )
            shutil.copy2(
                commandfile,
                commandfile + ts,
            )

        if not knownargs.nextflow:
            min_version("7.32.3")
            run_snakemake(
                knownargs.configfile,
                knownargs.directory,
                knownargs.use_conda,
                knownargs.procs,
                knownargs.skeleton,
                knownargs.loglevel,
                knownargs.save,
                knownargs.unlock,
                optionalargs[0],
            )

        else:
            nf_min_version = "23.04.1"
            nf_ver = mw.nf_check_version(nf_min_version)
            if nf_ver:
                run_nextflow(
                    knownargs.configfile,
                    knownargs.directory,
                    knownargs.use_conda,
                    knownargs.procs,
                    knownargs.skeleton,
                    knownargs.loglevel,
                    knownargs.save,
                    knownargs.clean,
                    optionalargs[0],
                )
            else:
                log.error(
                    logid
                    + "Minimal version of nextflow required is "
                    + str(nf_min_version)
                    + " and we only found "
                    + str(nf_ver)
                    + "! Please install or use envs/MONSDA.yaml to create conda environment accordingly"
                )
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


####################
####    MAIN    ####
####################


if __name__ == "__main__":
    logid = scriptname + ".main: "
    try:
        main()

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))

# MONSDA.py ends here
