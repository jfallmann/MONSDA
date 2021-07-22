#!/usr/bin/env python
from setuptools import setup, find_packages
from glob import glob
import os.path

NAME = "NextSnakes"
# Set __version__
exec(open('NextSnakes/__init__.py').read())
DESCRIPTION = "NextSnakes, a modular assembler of snakemake and nexflow workflows"

scripts = list()
for s in glob('scripts/*'):
    scripts.append(s)#os.path.join(s, os.path.split(s)[1]))
for d in glob('workflows/*'):
    scripts.append(d)#os.path.join(d, os.path.split(d)[1]))

requires = open("requirements.txt").read().strip().split("\n")


setup(
    name=NAME,
    version=__version__,
    description=DESCRIPTION,
    author="Joerg Fallmann",
    author_email="fall@bioinf.uni-leipzig.de",
    packages=find_packages(include=['NextSnakes', 'NextSnakes.*']),
    include_package_data=True,
    scripts=scripts,
    entry_points={
        "console_scripts": [
            "NextSnakes = NextSnakes.RunNextSnakes:main",
            "NextSnakes configure = NextSnakes.Configurator:main"
        ]
    },
    #scripts=["NextSnakes/Configurator.py"],
    install_requires=requires,
    python_requires=">=3.6",
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    zip_safe=False,
    data_files=[("", ["LICENSE"])],
    license='LICENSE',
    url="https://github.com/jfallmann/NextSnakes",
    long_description_content_type="text/markdown",
    long_description=open('README.md').read()
)
