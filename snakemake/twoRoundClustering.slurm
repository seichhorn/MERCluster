#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p zhuang,shared
#SBATCH -t 5000
#SBATCH --mem 1000
#SBATCH --open-mode=append
#SBATCH -o /n/home13/seichhorn/Scratch/%j.stdout.txt
#SBATCH -e /n/home13/seichhorn/Scratch/%j.stdout.txt


date +'Starting at %R.'

# snakemake command

/n/home13/seichhorn/.conda/envs/imagevirtualenv/bin/snakemake --jobs 400 -s /n/home13/seichhorn/Python_code/MERCluster/snakemake/twoRoundClustering.Snakefile --cluster-config /n/home13/seichhorn/Python_code/MERCluster/snakemake/twoRoundClustering_cluster.json --configfile /n/home13/seichhorn/Python_code/MERCluster/snakemake/twoRoundClustering_config.json --cluster "sbatch --mem={cluster.mem} -A {cluster.account} -p {cluster.partition} -n {cluster.n}  -t {cluster.time} -o {cluster.out} -e {cluster.err} --exclude={cluster.exclude} " /n/boslfs/LABS/zhuang_lab/User/seichhorn/Hypothalamus/SN-seq_180912/190107/Round2/analysis/cellTypes.txt

date +'Finished at %R.'
