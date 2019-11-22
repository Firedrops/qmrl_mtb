#!/bin/bash
#SBATCH --job-name=TBpipe
#SBATCH --output=TBpipe.%j.stdout
#SBATCH --error=TBpipe.%j.stderr
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32000MB
echo $SLURM_ARRAY_TASK_ID
if [ -z ${SLURM_ARRAY_TASK_ID+x} ]; then 
	echo "SLURM_ARRAY_TASK_ID is unset, setting it to 1";
        export SLURM_ARRAY_TASK_ID=1	
#	exit 1; 
else echo "SLURM_ARRAY_TASK_ID is set to '$SLURM_ARRAY_TASK_ID'"; 
fi
n=$SLURM_ARRAY_TASK_ID
pwd=$(pwd)
nme=$(ls ${pwd}/data | grep _R1.fastq.gz | sort | head -n $n  | tail -1 | cut -f 1 -d '_')
echo "running "$nme
##USE FOLLOWING COMMAND TO LOAD DOCKER IMAGE WITHOUT RUNNING SCRIPT FOR DEBUGGING
#docker run -it --entrypoint /bin/bash -v ${pwd}/data:/data/ dockersubtest/qimr_mtb:slim

#FOLLOWING COMMAND RUNS THE SCRIPT
docker run --rm -v ${pwd}/data:/data/ dockersubtest/qimr_mtb:slim /data/ $nme /data/H37Rv_refe /data/temp /data/out
