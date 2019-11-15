#!/bin/bash
#SBATCH --job-name=TBpipe
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16000 # mb
#SBATCH --time=9999:00:00

n=$SLURM_ARRAY_TASK_ID
nme=$(head -n $n /home/larry/png_mtb/filesin.txt | tail -1)

docker run -it --rm -v /home/larry/png_mtb/:/data/ dockersubtest/qimr_mtb:slim /data/ $nme /data/H37Rv_refe /data/temp /data/out
