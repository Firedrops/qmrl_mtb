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





#This command is to be run when all pairs have been processed, using data from all of them. WIP.

bcftools view -h ${tempdir}${NAME} | grep "^#CHROM" | cut -f10-
bcftools view -c1 -s ${NAME} -o ${tempdir}${NAME}_in.vcf ${tempdir}${NAME}.recode.vcf

python vcf_filter_module.py 9 ${tempdir}${NAME}_in.vcf ${tempdir}${NAME}_out.vcf

#Commands below are to be run on the file generated in the previous step.
#gatk3 -T CombineVariants -R /data/reference.fasta -–variant *.out.vcf -o master.vcf -genotypeMergeOptions UNIQUIFY
## shortcut alternative for when testing with 1 file:
#cat ${tempdir}${NAME}_out.vcf | bgzip -c > ${outdir}${NAME}_master.vcf.gz

#zcat ${outdir}${NAME}_master.vcf.gz | vcf-to-tab > ${tempdir}${NAME}_snps.tab

#/vcf_tab_to_fasta_alignment.pl -i ${tempdir}${NAME}_snps.tab > ${outdir}${NAME}_all_snps.fasta

#Make sure nectar users can access
chmod -R 777 ${tempdir}
chmod -R 777 ${outdir}

printf "Pipeline finished. Please check ${tempdir} and delete if not needed."

