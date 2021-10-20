VERSION=$1
#tag
git tag -f v$VERSION
#build
rm -rf .eggs build *.egg-info dist ; nocorrect python setup.py bdist_wheel sdist
#goto test dir
cd ~/Work/Test/Pipi
#uninstall old and install local new
pip uninstall -y NextSnakes && pip install /home/fall/NextSnakes/dist/NextSnakes-$VERSION\-py3-none-any.whl
#run snakemake
clear && NextSnakes -j 4 --configfile multitool.json --directory ${PWD} --conda-frontend mamba
#run nextflow
clear && NextSnakes --nextflow -j 4 -resume --configfile multitool.json --directory ${PWD}
