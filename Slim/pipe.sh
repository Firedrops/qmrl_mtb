#!/bin/bash

#Example usage: ./pipe.sh /data/ RB17MT1804 /data/H37Rv_refe /data/temp /data/out

#establishes variables
indir=$1
shift
NAME=$1
shift
reference=$1
shift
tempdir=$1_${NAME}/
shift
outdir=$1_${NAME}/
shift

#checks for pre-existing indices/dictionary
printf "Checking for bwa index.\nNote: Only checks for .amb and assumes other files are also present if found.\n"
if [ -f "${reference}.fasta.amb" ]
then
	printf "index found, skipping indexing step.\n"
else
	printf "index not found, proceeding with generating index.\n"
  printf "bwa index ${reference}.fasta"
  bwa index ${reference}.fasta
fi

printf "Checking for samtools index.\n"
if [ -f "${reference}.fasta.fai" ]
then
	printf "index found, skipping indexing step.\n"
else
	printf "index not found, proceeding with generating index.\n"
  printf "samtools faidx ${reference}.fasta"
  samtools faidx ${reference}.fasta
fi

printf "Checking for picard dictionary.\n"
if [ -f "${reference}.fasta.dict" ]
then
	printf "dictionary found, skipping generation step.\n"
else
	printf "dictionary not found, proceeding with generating dictionary.\n"
  printf "picard CreateSequenceDictionary R=${reference}.fasta O=${reference}.fasta.dict\n"
  picard CreateSequenceDictionary R=${reference}.fasta O=${reference}.fasta.dict
fi

printf "Copying *.fasta.dict to *.dict so that both are available.\n"
if [ -f "${reference}.dict" ]
then
	printf "*.dict found, skipping copying step.\n"
else
	printf "*.dict not found, proceeding with copying dictionary.\n"
  printf "picard CreateSequenceDictionary R=${reference}.fasta O=${reference}.fasta.dict\n"
  cp ${reference}.fasta.dict ${reference}.dict
fi

#make unique temporary and output directories
mkdir $tempdir
mkdir $outdir

#SKIPPED
#java -jar /trimmomatic/classes/trimmomatic.jar PE -phred33 -trimlog ${tempdir}${NAME}_log.txt ${indir}${NAME}_R1.fastq.gz ${indir}${NAME}_R2.fastq.gz ${tempdir}${NAME}_paired_R1.fastq.gz ${tempdir}${NAME}_unpaired_R1.fastq.gz ${tempdir}${NAME}_paired_R2.fastq.gz ${tempdir}${NAME}_unpaired_R2.fastq.gz ILLUMINACLIP:/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36
#SKIPPED

bwa mem -t 4 -M -R "@RG\tID:${NAME}\tSM:${NAME}\tPL:Illumina\tLB:001\tPU:001" ${reference}.fasta ${indir}${NAME}_R1.fastq.gz ${indir}${NAME}_R2.fastq.gz > ${tempdir}${NAME}.sam

samtools view -bS -F 4 ${tempdir}${NAME}.sam > ${tempdir}${NAME}.bam

picard SortSam I=${tempdir}${NAME}.bam O=${tempdir}${NAME}_sorted.bam SO=coordinate

picard MarkDuplicates I=${tempdir}${NAME}_sorted.bam O=${tempdir}${NAME}_dup.bam M=${tempdir}${NAME}.txt

picard BuildBamIndex I=${tempdir}${NAME}_dup.bam

gatk3 -T RealignerTargetCreator -R ${reference}.fasta -I ${tempdir}${NAME}_dup.bam -o ${tempdir}${NAME}.intervals

gatk3 -T IndelRealigner -R ${reference}.fasta -I ${tempdir}${NAME}_dup.bam -targetIntervals ${tempdir}${NAME}.intervals -o ${outdir}${NAME}_dup_alig.bam

#Make sure nectar users can access
chmod -R 777 ${tempdir}
chmod -R 777 ${outdir}

printf "Pipeline finished. Please check ${tempdir} and delete if not needed."
