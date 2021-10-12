#!/usr/bin/env python
from setuptools import setup, find_packages
from glob import glob
import os
import versioneer

NAME = "NextSnakes"
DESCRIPTION = "NextSnakes, a modular assembler of snakemake and nexflow workflows"
# Set __version__ done by versioneer
# exec(open("NextSnakes/__init__.py").read())

scripts = list()
for s in glob("scripts/**", recursive=True):
    if any(x in s for x in [".pl", ".py", ".sh"]):
        scripts.append(os.path.relpath(s))  # os.path.join(s, os.path.split(s)[1]))
workflows = list()
for d in glob("workflows/*"):
    if not "wip" in d:
        workflows.append(os.path.relpath(d))  # os.path.join(d, os.path.split(d)[1]))
envs = list()
for e in glob("envs/*"):
    envs.append(os.path.relpath(e))  # os.path.join(d, os.path.split(d)[1]))

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
    data_files=[
        (os.path.join("share", "NextSnakes", "workflows"), workflows),
        (os.path.join("share", "NextSnakes", "scripts"), scripts),
        (os.path.join("share", "NextSnakes", "envs"), envs),
        ("", ["LICENSE"]),
    ],
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
