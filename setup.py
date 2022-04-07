#!/usr/bin/env python

from setuptools import setup, find_packages
from glob import glob
from collections import defaultdict
import os
import versioneer


NAME = "MONSDA"
DESCRIPTION = (
    "MONSDA, Modular Organizer of Nextflow and Snakemake driven hts Data Analysis"
)
# Set __version__ done by versioneer
# exec(open("MONSDA/__init__.py").read())


def generate_datafiles():
    df = list()
    dirlist = defaultdict(list)

    scripts = list()
    for s in glob("scripts/**", recursive=True):
        if any(x in s for x in [".pl", ".pm", ".py", ".sh", ".R"]):
            scripts.append(os.path.relpath(s))  # os.path.join(s, os.path.split(s)[1]))
    for s in scripts:
        dirlist[str(os.path.join("share", "MONSDA", os.path.dirname(s)))].append(s)

    workflows = list()
    for d in glob("workflows/*"):
        if not "wip" in d:
            workflows.append(
                os.path.relpath(d)
            )  # os.path.join(d, os.path.split(d)[1]))
    for w in workflows:
        dirlist[os.path.join("share", "MONSDA", os.path.dirname(w))].append(w)

    envs = list()
    for e in glob("envs/*"):
        envs.append(os.path.relpath(e))  # os.path.join(d, os.path.split(d)[1]))
    for e in envs:
        dirlist[os.path.join("share", "MONSDA", os.path.dirname(e))].append(e)

    confs = list()
    for c in glob("configs/*"):
        if any(x in c for x in [".json"]):
            confs.append(os.path.relpath(c))  # os.path.join(d, os.path.split(d)[1]))
    for c in confs:
        dirlist[os.path.join("share", "MONSDA", os.path.dirname(c))].append(c)

    profiles = list()
    for p in glob("profile_*/*"):
        confs.append(os.path.relpath(p))  # os.path.join(d, os.path.split(d)[1]))
    for p in profiles:
        dirlist[os.path.join("share", "MONSDA", os.path.dirname(p))].append(p)

    dirlist[""].append("LICENSE")

    for k, v in dirlist.items():
        df.append((k, v))

    return df


# requires = open(os.path.abspath("requirements.txt")).read().strip().split("\n")

setup(
    name=NAME,
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description=DESCRIPTION,
    author="Joerg Fallmann",
    author_email="fall@bioinf.uni-leipzig.de",
    packages=find_packages(include=["MONSDA", "MONSDA.*"]),
    include_package_data=True,
    data_files=generate_datafiles(),
    entry_points={
        "console_scripts": [
            "monsda = MONSDA.RunMONSDA:main",
            "monsda_configure = MONSDA.Configurator:main",
        ]
    },
    # install_requires=requires,
    install_requires=[
        "biopython>=1.78",
        "snakemake>=6.5.3",
        "black>=21.5b2",
        "flake8>=3.8.3",
        "isort>=5.9.2",
        "sphinx>=4.1.0",
    ],
    python_requires=">=3.9",
    setup_requires=["pytest-runner"],
    tests_require=["pytest"],
    zip_safe=False,
    license="LICENSE",
    url="https://github.com/jfallmann/MONSDA",
    long_description_content_type="text/markdown",
    long_description=open("README.md").read(),
)
