# Versions.py ---

# Export of tool versions of a MONSDA run as deduced from the config file and
# the conda environment yaml files of the configured tools.

import datetime
import logging
import os
import re

log = logging.getLogger(__name__)

scriptname = os.path.basename(__file__).replace(".py", "")

# Tools/packages that are no analysis tools and only clutter the version list
SKIP_PACKAGES = frozenset(
    [
        "pip",
        "python",
        "setuptools",
        "wheel",
        "openssl",
        "zlib",
    ]
)

# Mapping of tool/package name to the corresponding MultiQC module name(s),
# names have to match the modules of MultiQC exactly, including case,
# tools without a MultiQC module are not listed, tools reporting the output of
# several other tools, like rustqc, map to all modules they produce output for
MULTIQC_MODULES = {
    "bbmap": "bbmap",
    "bbduk": "bbduk",
    "bismark": "bismark",
    "bowtie2": "bowtie2",
    "cutadapt": "cutadapt",
    "fastp": "fastp",
    "fastqc": "fastqc",
    "featurecounts": "featureCounts",
    "hisat2": "hisat2",
    "htseq": "htseq",
    "kallisto": "kallisto",
    "macs2": "macs2",
    "macs3": "macs2",
    "picard": "picard",
    "qualimap": "qualimap",
    "rastqc": "fastqc",
    "rseqc": "rseqc",
    "rustqc": (
        "custom_content",
        "featureCounts",
        "preseq",
        "qualimap",
        "rseqc",
        "samtools",
    ),
    "salmon": "salmon",
    "samtools": "samtools",
    "star": "star",
    "subread": "featureCounts",
    "trimgalore": "cutadapt",
    "trim-galore": "cutadapt",
    "trimmomatic": "trimmomatic",
    "umitools": "umitools",
    "umi_tools": "umitools",
}


def multiqc_module(tool):
    """Return the MultiQC module(s) for a tool or None if there is none

    Parameters
    ----------
    tool : str
        name of tool/conda package

    Returns
    -------
    str or None
        module name, comma separated if the tool reports for several modules
    """
    module = MULTIQC_MODULES.get(
        str(tool).lower().replace("_", "").replace("-", "")
    ) or MULTIQC_MODULES.get(str(tool).lower())
    if isinstance(module, tuple):
        return ",".join(module)
    return module


def envs_from_config(config, steps):
    """Collect conda environments of all tools configured to run

    Parameters
    ----------
    config : dict
        config json file parsed by snakemake.common.configfile
    steps : list
        list of workflow steps (e.g. ['QC', 'TRIMMING', 'MAPPING']) that will be run

    Returns
    -------
    dict
        environment name to set of workflow steps it is used in
    """
    logid = scriptname + ".envs_from_config: "
    envs = dict()

    for step in steps:
        if not step or not config.get(step):
            continue
        conf = config[step]
        if not isinstance(conf, dict):
            continue
        found = set()
        stack = [conf]
        while stack:
            sub = stack.pop()
            for key, val in sub.items():
                if isinstance(val, dict):
                    if key == "TOOLS":
                        found.update(
                            str(x) for x in val.keys() if x and str(x) != "None"
                        )
                    else:
                        stack.append(val)
                elif key.endswith("ENV") and isinstance(val, str) and val:
                    found.add(str(val))
        for env in found:
            envs.setdefault(env, set()).add(step)

    log.debug(logid + str(envs))
    return envs


def env_yaml(env, step, envdir):
    """Find the conda env yaml belonging to an env name of a workflow step

    Parameters
    ----------
    env : str
        environment name as configured
    step : str
        workflow step the env is used in
    envdir : str
        directory holding the conda env yamls

    Returns
    -------
    str or None
        path to yaml file
    """
    for name in [env, "_".join([env, step]), "_".join([env, str(step).lower()])]:
        path = os.path.join(envdir, name + ".yaml")
        if os.path.isfile(path):
            return path
    return None


def parse_env_yaml(yamlfile):
    """Parse pinned versions from the dependencies of a conda env yaml

    Parameters
    ----------
    yamlfile : str
        path to conda env yaml

    Returns
    -------
    dict
        package name to version string, 'unpinned' if no version is defined
    """
    logid = scriptname + ".parse_env_yaml: "
    versions = dict()
    indeps = False

    with open(yamlfile, "r") as y:
        for line in y:
            line = line.rstrip()
            if not line or line.lstrip().startswith("#"):
                continue
            entry = line.strip()
            if not entry.startswith("- "):
                indeps = line.startswith("dependencies:")
                continue
            if not indeps:
                continue
            entry = entry[2:].strip()
            if entry.endswith(":"):
                continue
            match = re.match(r"^([A-Za-z0-9._+-]+)\s*(.*)$", entry)
            if not match:
                continue
            pkg, spec = match.group(1), match.group(2).strip()
            if pkg.lower() in SKIP_PACKAGES:
                continue
            spec = re.sub(r"^[=]{1,2}\s*", "", spec)
            versions[pkg] = spec if spec else "unpinned"

    log.debug(logid + yamlfile + ": " + str(versions))
    return versions


# Package name prefixes of conda channels that are no part of the tool name
PACKAGE_PREFIXES = ("bioconductor-", "r-", "python-", "perl-")

# Tools always reported if part of an env, independent of the env name
ALWAYS_REPORT = frozenset(["multiqc"])


def _norm(name):
    name = str(name).lower()
    for prefix in PACKAGE_PREFIXES:
        if name.startswith(prefix) and len(name) > len(prefix):
            name = name[len(prefix) :]
            break
    return re.sub(r"[_\-.]", "", name)


def relevant_packages(env, packages):
    """Reduce the packages of an env to the tools worth reporting

    Keeps the package the env is named after, plus every package MultiQC has a
    module for, so dependencies like perl or matplotlib do not clutter the list.

    Parameters
    ----------
    env : str
        environment name as configured
    packages : dict
        package name to version as parsed from the env yaml

    Returns
    -------
    dict
        package name to version
    """
    keep = dict()
    nenv = _norm(env)

    for pkg, version in packages.items():
        npkg = _norm(pkg)
        if npkg == nenv or npkg.startswith(nenv) or nenv.startswith(npkg):
            keep[pkg] = version
        elif multiqc_module(pkg) or pkg.lower() in ALWAYS_REPORT:
            keep[pkg] = version

    return keep


def collect_versions(config, steps, envdir):
    """Collect versions of all tools configured to run

    Parameters
    ----------
    config : dict
        config json file parsed by snakemake.common.configfile
    steps : list
        list of workflow steps that will be run
    envdir : str
        directory holding the conda env yamls

    Returns
    -------
    list
        rows of [tool, version, multiqc module, steps] sorted by tool
    """
    logid = scriptname + ".collect_versions: "
    tools = dict()

    for env, envsteps in envs_from_config(config, steps).items():
        yamlfile = None
        for step in sorted(envsteps):
            yamlfile = env_yaml(env, step, envdir)
            if yamlfile:
                break
        if not yamlfile:
            log.warning(
                logid + "No conda env yaml found for " + str(env) + ", skipping version"
            )
            tools.setdefault(env, ["unknown", set()])[1].update(envsteps)
            continue
        packages = relevant_packages(env, parse_env_yaml(yamlfile))
        if not packages:
            log.warning(
                logid
                + "No version for "
                + str(env)
                + " deducible from "
                + str(yamlfile)
            )
            packages = {env: "unknown"}
        # tools shipped by a differently named package, e.g. bbduk by bbmap,
        # need their own MultiQC module to be reported
        envmod = multiqc_module(env)
        if envmod and envmod not in [multiqc_module(pkg) for pkg in packages]:
            packages[env] = (
                list(packages.values())[0] if len(packages) == 1 else "unknown"
            )
        for pkg, version in packages.items():
            tools.setdefault(pkg, [version, set()])
            tools[pkg][0] = version
            tools[pkg][1].update(envsteps)

    rows = list()
    for tool in sorted(tools):
        version, envsteps = tools[tool]
        rows.append(
            [
                tool,
                version,
                multiqc_module(tool) or "-",
                ",".join(sorted(envsteps)),
            ]
        )

    log.debug(logid + str(rows))
    return rows


def write_versions(config, steps, envdir, monsda_version, configfile=None, path=None):
    """Write versions of all configured tools to LOGS/versions.txt

    Parameters
    ----------
    config : dict
        config json file parsed by snakemake.common.configfile
    steps : list
        list of workflow steps that will be run
    envdir : str
        directory holding the conda env yamls
    monsda_version : str
        version of MONSDA
    configfile : str, optional
        config file name to note in the header
    path : str, optional
        output file, by default LOGS/versions.txt

    Returns
    -------
    str
        path of the versions file
    """
    logid = scriptname + ".write_versions: "
    path = path if path else os.path.join("LOGS", "versions.txt")
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)

    if not os.path.isdir(envdir):
        fallback = os.path.join(os.path.dirname(os.path.dirname(__file__)), "envs")
        log.warning(
            logid + str(envdir) + " not found, falling back to " + str(fallback)
        )
        envdir = fallback

    rows = collect_versions(config, steps, envdir)

    with open(path, "w") as out:
        out.write(
            "# MONSDA {v}\t{ts}\n".format(
                v=monsda_version,
                ts=datetime.datetime.now().strftime("%Y%m%d_%H_%M_%S"),
            )
        )
        if configfile:
            out.write("# config: {c}\n".format(c=str(configfile)))
        out.write(
            "# versions as pinned in the conda env yamls of the configured tools, not queried from the running binaries\n"
        )
        out.write("#tool\tversion\tmultiqc_module\tsteps\n")
        for row in rows:
            out.write("\t".join(row) + "\n")

    log.info(logid + "Tool versions written to " + str(path))
    return path
