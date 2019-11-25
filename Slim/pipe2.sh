#!/bin/bash
#Example usage: ./pipe.sh /data/ /data/H37Rv_refe /data/temp /data/out

#establishes variables
DATE=$(date +"%Y%m%d")
NUM=$(ls out/out_*/*_dup_alig.bam | wc -l )
NAME=${NME}_${DATE}_${NUM} (edited) 
input=$(ls out/out_*/*_dup_alig.bam)

indir=$1
shift
reference=$1
shift
tempdir=$1_${NAME}/
shift
outdir=$1_${NAME}/
shift



gatk3 -T UnifiedGenotyper -R ${reference}.fasta -I ${input} -A AlleleBalance -pnrm EXACT_GENERAL_PLOIDY -ploidy 1 -glm SNP -o ${outdir}${NAME}.vcf

gatk3 -T VariantFiltration -R ${reference}.fasta -V ${outdir}${NAME}.vcf --filterExpression "((DP-MQ0)<10) || ((MQ0/(1.0*DP))>=0.8) || (ABHom <0.8) || (Dels >0.5) || (QUAL > 90)" --filterName LowConfidence -o ${outdir}${NAME}_filtered.vcf

vcftools --vcf ${outdir}${NAME}_filtered.vcf --out ${tempdir}${NAME} --recode --keep-INFO-all

python vcf_filter_module.py 9 ${tempdir}${NAME}_in.vcf ${tempdir}${NAME}_out.vcf

gatk3 -T CombineVariants -R ${reference}.fasta -–variant ${tempdir}${NAME}_out.vcf -o ${outdir}${NAME}_master.vcf -genotypeMergeOptions UNIQUIFY


zcat ${outdir}${NAME}_master.vcf | vcf-to-tab > ${tempdir}${NAME}_snps.tab

perl /vcf_tab_to_fasta_alignment.pl -i ${tempdir}${NAME}_snps.tab > ${outdir}${NAME}_all_snps.fasta

#Make sure nectar users can access
chmod -R 777 ${tempdir}
chmod -R 777 ${outdir}

printf "Pipeline finished. Please check ${tempdir} and delete if not needed."

