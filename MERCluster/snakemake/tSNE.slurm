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




date +'Finished at %R.'
