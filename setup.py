from setuptools import setup

setup(
    name='NextSnakes',
    version='3.0.0',
    author='Joerg Fallmann',
    author_email='fall@bioinf.uni-leipzig.de',
    packages=find_packages(),
    scripts=[],
    license='LICENSE',
    url="https://github.com/jfallmann/NextSnakes",
    long_description_content_type="text/markdown",
    description='',
    long_description=open('README.md').read(),
    include_package_data=True,
    install_requires=[
        "biopython>=1.78",
        "snakemake>=6.5.3",
        "black>=21.5",
        "flake8>=3.8.3",
        "isort>=5.5.2",
        "sphinx>=4.1.0",
    ],
    python_requires=">=3.6",
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
)
