#!/usr/bin/env python
from setuptools import setup, find_packages
from glob import glob
import os
import versioneer

NAME = "NextSnakes"
DESCRIPTION = "NextSnakes, a modular assembler of snakemake and nexflow workflows"
# Set __version__ done by versioneer
# exec(open("NextSnakes/__init__.py").read())

def generate_datafiles():
    df = list()
    scripts = list()
    for s in glob("scripts/**", recursive=True):
        if any(x in s for x in [".pl", ".py", ".sh"]):
            scripts.append(os.path.relpath(s))  # os.path.join(s, os.path.split(s)[1]))
    for s in scripts:
        df.append((os.path.join("share", "NextSnakes", os.path.dirname(s)), [os.path.basename(s)]))

    workflows = list()
    for d in glob("workflows/*"):
        if not "wip" in d:
            workflows.append(os.path.relpath(d))  # os.path.join(d, os.path.split(d)[1]))
    for w in workflows:
        df.append((os.path.join("share", "NextSnakes", os.path.dirname(w)), [os.path.basename(w)]))

    envs = list()
    for e in glob("envs/*"):
        envs.append(os.path.relpath(e))  # os.path.join(d, os.path.split(d)[1]))
    for e in envs:
        df.append((os.path.join("share", "NextSnakes", os.path.dirname(e)), [os.path.basename(e)]))

    df.append(("", ["LICENSE"]))
    print(df)
    return(df)

requires = open("requirements.txt").read().strip().split("\n")

setup(
    name=NAME,
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description=DESCRIPTION,
    author="Joerg Fallmann",
    author_email="fall@bioinf.uni-leipzig.de",
    packages=find_packages(include=["NextSnakes", "NextSnakes.*"]),
    include_package_data=True,
    # scripts=scripts,
    data_files=generate_datafiles(),
    entry_points={
        "console_scripts": [
            "NextSnakes = NextSnakes.RunNextSnakes:main",
            "NextSnakes_configure = NextSnakes.Configurator:main",
        ]
    },
    # scripts=["NextSnakes/Configurator.py"],
    install_requires=requires,
    python_requires=">=3.6",
    setup_requires=["pytest-runner"],
    tests_require=["pytest"],
    zip_safe=False,
    license="LICENSE",
    url="https://github.com/jfallmann/NextSnakes",
    long_description_content_type="text/markdown",
    long_description=open("README.md").read(),
)
