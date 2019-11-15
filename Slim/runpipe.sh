#!/bin/bash
#SBATCH --job-name=TBpipe
#SBATCH --output=TBpipe.%j.stdout
#SBATCH --error=TBpipe.%j.stderr
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32000MB

n=$SLURM_ARRAY_TASK_ID
nme=$(head -n $n filesin.txt | tail -1 | cut -f 1 -d '_')

docker run -it --rm -v ${pwd}/:/data/ dockersubtest/qimr_mtb:slim /data/ $nme /data/H37Rv_refe /data/temp /data/out
