from setuptools import setup, find_packages

NAME = "NextSnakes"
VERSION = "1.0.6"
DESCRIPTION = "NextSnakes, a modular assembler of snakemake and nexflow workflows"

setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    author="Joerg Fallmann",
    author_email="fall@bioinf.uni-leipzig.de",
    packages=find_packages(),
    scripts=["NextSnakes.py", "Configurator.py"],
    license='LICENSE',
    url="https://github.com/jfallmann/NextSnakes",
    long_description_content_type="text/markdown",
    long_description=open('README.md').read(),
    include_package_data=True,
    install_requires=[
        "biopython>=1.78",
        "snakemake>=6.5.3",
        "black>=21.5b2",
        "flake8>=3.8.3",
        "isort>=5.9.2",
        "sphinx>=4.1.0",
    ],
    python_requires=">=3.6",
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
)
