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
##IF DEBUGGING
#NAME=20191201_96; indir="data/"; reference="H37Rv"; tempdir=${indir}_${NAME}/ ; outdir=${indir}_${NAME}/
tempdir=${indir}_${NAME}/
#shift
outdir=${indir}_${NAME}/
#shift




mkdir ${outdir}

gatk3 -T UnifiedGenotyper -R ${reference}.fasta -I /out/bams.list -A AlleleBalance -pnrm EXACT_GENERAL_PLOIDY -ploidy 1 -glm SNP -o ${outdir}${NAME}.vcf

gatk3 -T VariantFiltration -R ${reference}.fasta -V ${outdir}${NAME}.vcf --filterExpression "((DP-MQ0)<10) || ((MQ0/(1.0*DP))>=0.8) || (ABHom <0.8) || (Dels >0.5) || (QUAL > 90)" --filterName LowConfidence -o ${outdir}${NAME}_filtered.vcf

pwd=$(pwd)
cd $outdir

## this outputs as out_recode.vcf
 vcftools --vcf ${NAME}_filtered.vcf  --recode --keep-INFO-all --remove-filtered-all

#Currently redundant
##FOLLOWING COMMAND SHOULD 

file="out.recode.vcf"
for sample in $(bcftools view -h $file | grep "^#CHROM" | cut -f10-); do  bcftools view -c1 -s $sample -o ${file/.vcf*/.$sample.vcf} $file; done

for file1 in $(ls out.recode.*vcf | grep -v 'master.vcf'  | grep -v out.recode.vcf); do

#python vcf_filter_module.py 9 ${tempdir}${NAME}_filtered.vcf ${tempdir}${NAME}_master.vcf
outfile1=${file1}_master.vcf
python ${pwd}/vcf_filter_module.py 12 ${file1} ${outfile1}
#Currently redundant. Note: Rename ${tempdir}${NAME}_master.vcf to ${tempdir}${NAME}_out.vcf above in line 34 if this command is re-instated.
done


 
 a=$(ls out.recode.*_master.vcf |  sed 's/out.recode/--variant out.recode/g')

gatk3 -T CombineVariants -R ${pwd}/${reference}.fasta  ${a} -o ${NAME}_master2.vcf -genotypeMergeOptions UNIQUIFY


cat ${NAME}_master2.vcf | vcf-to-tab > ${NAME}_snps.tab


##THIS IS ATTEMPTING TO RUN SNPEFF
##NOTE THE SNPEFF COMMAND DOES NOT SEEM TO BE WORKING
cat ${NAME}_master2.vcf | sed -e 's/NC_000962.3/NC_000962/' >${NAME}_master3.vcf
java -jar  /snpEff/snpEff.jar -c ../snpEff.config -v ../m_tuberculosis_H37Rv > ${NAME}_master3.vcf > ${NAME}_master4.vcf



#FOLLOWING LINE REPLACES DOT WITH REFERENCE FROM THAT POSITION
while read line; do ref=$(echo $line | cut -f 3 -d ' '); echo $line | sed "s/\.variant//g" | sed "s/\./$ref/g"; done <  ${NAME}_snps.tab  > ${NAME}_snps1.tab
sed 's/ /\t/g' ${NAME}_snps1.tab  > ${NAME}_snps2.tab
 
perl ${pwd}/vcf_tab_to_fasta_alignment.pl -i ${NAME}_snps2.tab > ${NAME}_all_snps.fasta
cd $pwd

#Make sure nectar users can access
chmod -R 777 ${tempdir}
chmod -R 777 ${outdir}

printf "Pipeline finished. Please check ${tempdir} and delete if not needed."

