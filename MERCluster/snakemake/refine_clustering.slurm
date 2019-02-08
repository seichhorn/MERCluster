#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p zhuang,shared
#SBATCH -t 5000
#SBATCH --mem 100000
#SBATCH --open-mode=append
#SBATCH -o /n/home13/seichhorn/Scratch/%j.stdout.txt
#SBATCH -e /n/home13/seichhorn/Scratch/%j.stdout.txt


/n/home13/seichhorn/.conda/envs/imagevirtualenv/bin/python /n/home13/seichhorn/Python_code/MERCluster/refining_clustering.py Astrocyte