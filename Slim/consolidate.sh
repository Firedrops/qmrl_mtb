  
#!/bin/bash

#Example usage: ./consolidate.sh /data/ RB17MT1804 /data/H37Rv_refe /data/temp /data/out

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

for dir in /data/* ; do
        realpath=$(ls $dir | grep filtered.vcf$) > /data/vcfin.txt;
done

for ${indir}${NAME} in ${indir}${NAME}s_filtered.vcf; 
  do for sample in `bcftools view -h ${indir}${NAME} | grep "^#CHROM" | cut  -f10-`; 
  do bcftools view -c1 -s $sample -o ${indir}${NAME}/.vcf*/.$sample.vcf} ${indir}${NAME};
  done;
done
  
## shortcut alternative for when testing with 1 file:
#bcftools view -c1 -s ${NAME} -o ${tempdir}${NAME}_in.vcf ${tempdir}${NAME}.recode.vcf

#Commands below are to be run on the file generated in the previous step.
#python vcf_filter_module.py 9 ${tempdir}${NAME}_in.vcf ${tempdir}${NAME}_out.vcf

#gatk3 -T CombineVariants -R /data/reference.fasta -–variant *.out.vcf -o master.vcf -genotypeMergeOptions UNIQUIFY
## shortcut alternative for when testing with 1 file:
#cat ${tempdir}${NAME}_out.vcf | bgzip -c > ${outdir}${NAME}_master.vcf.gz

#zcat ${outdir}${NAME}_master.vcf.gz | vcf-to-tab > ${tempdir}${NAME}_snps.tab

#/vcf_tab_to_fasta_alignment.pl -i ${tempdir}${NAME}_snps.tab > ${outdir}${NAME}_all_snps.fasta
