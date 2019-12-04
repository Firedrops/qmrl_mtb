#!/bin/bash
#Example usage: docker run -it --rm --entrypoint /pipe2.sh -v /home/lcoin/nextflow/data:/data/ -v /home/lcoin/nextflow/out:/out/ qimr_slim /data/ /data/H37Rv_refe /out/temp /out/out

#establishes variables
find /out -name *_dup_alig.bam > /out/bams.list
DATE=$(date +"%Y%m%d")
NUM=$(wc -l /out/bams.list  | cut -f 1 -d ' ')
NAME=${DATE}_${NUM}
#input=$(ls /out/out_*/*_dup_alig.bam)

indir=$1
shift
reference=$1
shift
tempdir=${indir}_${NAME}/
shift
outdir=${indir}_${NAME}/
shift




mkdir ${outdir}

gatk3 -T UnifiedGenotyper -R ${reference}.fasta -I /out/bams.list -A AlleleBalance -pnrm EXACT_GENERAL_PLOIDY -ploidy 1 -glm SNP -o ${outdir}${NAME}.vcf

gatk3 -T VariantFiltration -R ${reference}.fasta -V ${outdir}${NAME}.vcf --filterExpression "((DP-MQ0)<10) || ((MQ0/(1.0*DP))>=0.8) || (ABHom <0.8) || (Dels >0.5) || (QUAL > 90)" --filterName LowConfidence -o ${outdir}${NAME}_filtered.vcf

#Currently redundant
#bcftools view -h ${outdir}${NAME}_filtered.vcf | grep "^#CHROM" | cut -f10-

vcftools --vcf ${outdir}${NAME}_filtered.vcf --out ${tempdir}${NAME} --recode --keep-INFO-all

python vcf_filter_module.py 9 ${tempdir}${NAME}_filtered.vcf ${tempdir}${NAME}_master.vcf

#Currently redundant. Note: Rename ${tempdir}${NAME}_master.vcf to ${tempdir}${NAME}_out.vcf above in line 34 if this command is re-instated.
#gatk3 -T CombineVariants -R ${reference}.fasta -–variant ${tempdir}${NAME}_out.vcf -o ${outdir}${NAME}_master.vcf -genotypeMergeOptions UNIQUIFY

zcat ${outdir}${NAME}_master.vcf | vcf-to-tab > ${tempdir}${NAME}_snps.tab

perl /vcf_tab_to_fasta_alignment.pl -i ${tempdir}${NAME}_snps.tab > ${outdir}${NAME}_all_snps.fasta

#Make sure nectar users can access
chmod -R 777 ${tempdir}
chmod -R 777 ${outdir}

printf "Pipeline finished. Please check ${tempdir} and delete if not needed."

