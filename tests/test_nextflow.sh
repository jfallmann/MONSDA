DEF="INFO"
LVL="${1:-$DEF}"
bash cleanup.sh && export NXF_EXECUTOR=slurm; MONSDA --nextflow -j 8 --configfile multitool.json --directory ${PWD} --loglevel $LVL
