#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p zhuang,shared,unrestricted
#SBATCH -t 5000
#SBATCH --mem 1000
#SBATCH --open-mode=append
#SBATCH -o /n/home13/seichhorn/Scratch/%j.stdout.txt
#SBATCH -e /n/home13/seichhorn/Scratch/%j.stdout.txt


date +'Starting at %R.'

# snakemake command

/n/home13/seichhorn/.conda/envs/imagevirtualenv/bin/snakemake --jobs 400 -s /n/home13/seichhorn/Python_code/MERCluster/MERCluster/snakemake/singleRoundClustering.Snakefile --cluster-config /n/home13/seichhorn/Python_code/MERCluster/MERCluster/snakemake/twoRoundClustering_cluster.json --configfile /n/boslfs/LABS/zhuang_lab/User/seichhorn/PAG/merfish/clustering_200129/round1.json --cluster "sbatch --mem={cluster.mem} -A {cluster.account} -p {cluster.partition} -n {cluster.n}  -t {cluster.time} -o {cluster.out} -e {cluster.err} "

date +'Finished at %R.'
