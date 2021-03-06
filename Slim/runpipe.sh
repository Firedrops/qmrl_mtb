#!/bin/bash
#SBATCH --job-name=TBpipe
#SBATCH --output=TBpipe.%j.stdout
#SBATCH --error=TBpipe.%j.stderr
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32000MB
DEBUG=0

if [ -z ${SLURM_ARRAY_TASK_ID+x} ]; then 
	echo "SLURM_ARRAY_TASK_ID is unset, setting it to 1";
        export SLURM_ARRAY_TASK_ID=1	
#	exit 1; 
else 
	echo "SLURM_ARRAY_TASK_ID is set to '$SLURM_ARRAY_TASK_ID'"; 
fi
n=$SLURM_ARRAY_TASK_ID
pwd=$(pwd)
mkdir -p ${pwd}/out
dataline=$(readlink -f $pwd/data)
outputline=$(readlink -f $pwd/out)
nme=$(ls ${pwd}/data | grep _R1.fastq.gz | sort | head -n $n  | tail -1 | cut -f 1 -d '_')
echo "running "$nme

#docker images | grep slim | tr -s ' '  | cut -f 3 -d ' '  | xargs docker image rm
#THIS ASSUMES THERE IS A DOCKER SUBDIRECTOR TO BUILD THE DOCKER (IT IS A SYMBOLIC LINK TO Slim/ subdirectory
#cd docker/
#docker build -t qmrl_slim .
#cd $pwd
#docker pull dockersubtest/qmrl_mtb:slim

if [[ $DEBUG == 1 ]]; then
	##USE FOLLOWING COMMAND TO LOAD DOCKER IMAGE WITHOUT RUNNING SCRIPT FOR DEBUGGING
	docker run -it --entrypoint /bin/bash -v ${dataline}:/data/  -v ${outputline}:/out/  dockersubtest/qmrl_mtb:slim
else
	#FOLLOWING COMMAND RUNS THE SCRIPT
	#docker run --rm -v ${pwd}/:/data/ dockersubtest/qmrl_mtb:slim /data/ $nme /data/H37Rv_refe /data/temp /data/out
        docker run --rm -v ${dataline}:/data/ -v ${outputline}:/out/ dockersubtest/qmrl_mtb:slim /data/ $nme /data/H37Rv_refe /out/temp /out/out	

fi
