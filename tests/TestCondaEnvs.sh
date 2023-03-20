for i in ~/MONSDA/envs/*.yaml;do rm -rf ~/anaconda3/envs/tempenv;echo "INSTALLING $i" && conda env create -n tempenv -f $i && echo "DONE, NEXT";done
