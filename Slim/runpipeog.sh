#SBATCH --cpus-per-task=8
#SBATCH --mem=32000MB	#SBATCH --mem=32000MB

n=$SLURM_ARRAY_TASK_ID
pwd=$(pwd)
nme=$(ls ${pwd}/data | grep _R1.fastq.gz | sort | head -n $n  | tail -1 | cut -f 1 -d '_')
echo "running "$nme
docker run --rm -v ${pwd}/:/data/ dockersubtest/qimr_mtb:slim /data/ $nme /data/H37Rv_refe /data/temp /data/out
