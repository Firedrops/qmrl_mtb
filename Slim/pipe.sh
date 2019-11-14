#Example usage: ./pipe.sh /data/ RB17MT1804 /data/H37Rv_refe.fasta /data/temp /data/out

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

bwa index ${reference}
samtools faidx ${reference}
picard CreateSequenceDictionary R=${reference} O=${reference}.dict
mkdir $tempdir
mkdir $outdir

#java -jar /trimmomatic/classes/trimmomatic.jar PE -phred33 -trimlog ${tempdir}${NAME}_log.txt ${indir}${NAME}_R1.fastq.gz ${indir}${NAME}_R2.fastq.gz ${tempdir}${NAME}_paired_R1.fastq.gz ${tempdir}${NAME}_unpaired_R1.fastq.gz ${tempdir}${NAME}_paired_R2.fastq.gz ${tempdir}${NAME}_unpaired_R2.fastq.gz ILLUMINACLIP:/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36

bwa mem -t 4 -M -R "@RG\tID:${NAME}\tSM:${NAME}\tPL:Illumina\tLB:001\tPU:001" ${reference} ${indir}${NAME}_R1.fastq.gz ${indir}${NAME}_R2.fastq.gz > ${tempdir}${NAME}.sam

samtools view -bS -F 4 ${tempdir}${NAME}.sam > ${tempdir}${NAME}.bam

picard SortSam I=${tempdir}${NAME}.bam O=${tempdir}${NAME}_sorted.bam SO=coordinate

picard MarkDuplicates I=${tempdir}${NAME}_sorted.bam O=${tempdir}${NAME}_dup.bam M=${tempdir}${NAME}.txt

picard BuildBamIndex I=${tempdir}${NAME}_dup.bam

gatk3 -T RealignerTargetCreator -R ${reference} -I ${tempdir}${NAME}_dup.bam -o ${tempdir}${NAME}.intervals

gatk3 -T IndelRealigner -R ${reference} -I ${tempdir}${NAME}_dup.bam -targetIntervals ${tempdir}${NAME}.intervals -o ${outdir}${NAME}_dup_alig.bam

gatk3 -T UnifiedGenotyper -R ${reference} -I ${outdir}${NAME}_dup_alig.bam -A AlleleBalance -pnrm EXACT_GENERAL_PLOIDY -ploidy 1 -glm SNP -o ${outdir}${NAME}.vcf

gatk3 -T VariantFiltration -R ${reference} -V ${outdir}${NAME}.vcf --filterExpression "((DP-MQ0)<10) || ((MQ0/(1.0*DP))>=0.8) || (ABHom <0.8) || (Dels >0.5) || (QUAL > 90)" --filterName LowConfidence -o ${outdir}${NAME}_filtered.vcf

vcftools --vcf ${outdir}${NAME}_filtered.vcf --recode --keep-INFO-all

#for ${indir}${NAME} in ${indir}${NAME}s_filtered.vcf;do for sample in `bcftools view -h $${indir}${NAME} | grep "^#CHROM" | cut  -f10-`; do bcftools view -c1 -s $sample -o ${${indir}${NAME}/.vcf*/.$sample.vcf} $${indir}${NAME};done;done
bcftools view -c1 -s RB17MT0722 -o ${tempdir}${NAME}_in.vcf ${tempdir}${NAME}_out.recode.vcf

python vcf_filter_module.py 9 ${tempdir}${NAME}_in.vcf ${tempdir}${NAME}_out.vcf

#gatk3 -T CombineVariants -R /data/reference.fasta –V vcf1 vcf2 vcf3 -genotypeMergeOptions UNIQUIFY –o master.vcf
#gatk3 -T CombineVariants -R /data/reference.fasta -–variant out.vcf -o master.vcf -genotypeMergeOptions UNIQUIFY
cat ${tempdir}${NAME}_out.vcf | bgzip -c > ${outdir}${NAME}_master.vcf.gz

zcat ${outdir}${NAME}_master.vcf.gz | vcf-to-tab > ${tempdir}${NAME}_snps.tab

./vcf_tab_to_fasta_alignment.pl -i ${tempdir}${NAME}_snps.tab > ${outdir}${NAME}_all_snps.fasta
