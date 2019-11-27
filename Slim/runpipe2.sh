#!/bin/bash
#SBATCH --job-name=TBpipe2
#SBATCH --output=TBpipe2.%j.stdout
#SBATCH --error=TBpipe2.%j.stderr
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32000MB
DEBUG=0



pwd=$(pwd)
dataline=$(readlink -f $pwd/data)
outputline=$(readlink -f $pwd/out)
chmod a+x ${pwd}/pipe2.sh
slim=$(readlink -f $pwd/pipe2.sh  | rev | cut -f 2- -d / | rev)


if [[ $DEBUG == 1 ]]; then
	##USE FOLLOWING COMMAND TO LOAD DOCKER IMAGE WITHOUT RUNNING SCRIPT FOR DEBUGGING
	docker run -it --entrypoint /bin/bash -v ${dataline}:/data/  -v ${outputline}:/out/ dockersubtest/qimr_mtb:slim
else
	#FOLLOWING COMMAND RUNS THE SCRIPT
        docker run --rm --entrypoint /slim/pipe2.sh -v ${slim}:/slim/ -v ${dataline}:/data/ -v ${outputline}:/out/ dockersubtest/qimr_mtb:slim /data/ $nme /data/H37Rv_refe /out/temp /out/out	
fi
