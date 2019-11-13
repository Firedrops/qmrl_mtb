
java -jar /home/trimmomatic/current/bin/trimmomatic-0.32.jar PE -phred33 -trimlog ${file}_log.txt ${file}_R1.fastq.gz ${file}_R2.fastq.gz ${file}_paired_R1.fastq.gz ${file}_unpaired_R1.fastq.gz ${file}_paired_R2.fastq.g    z ${file}_unpaired_R2.fastq.gz ILLUMINACLIP:/home/trimmomatic/current/bin/adapters/NexteraPE-PE.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36


bwa mem -t 4 -M -R "@RG\tID:OOOO\tSM:${NAME}\tPL:Illumina\tLB:001\tPU:001" reference.fasta file_R1.fastq.gz file_R2.fastq.gz > file.sam

samtools view -bS -F 4 file.sam > file.bam

java -Xms1000M -Djava.io.tmpdir=./tmp. -jar picard-tools-1.119/SortSam.jar I=file.bam O=file_sorted.bam SO=coordinate


java -Xms100M -Djava.io.tmpdir=./tmp. -jar picard-tools-1.119/MarkDuplicates.jar I=file_sorted.bam O=file_dup.bam M=file.txt


java -Xms100M -Djava.io.tmpdir=./tmp. -jar picard-tools-1.119/BuildBamIndex.jar I=file_dup.bam


java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R reference.fasta -I file_dup.bam -o file.intervals


java -jar GenomeAnalysisTK.jar -T IndelRealigner -R reference.fasta -I file_dup.bam -targetIntervals file.intervals -o file_dup_alig.bam


java -jar GenomeAnalysisTK.jar -T UnifiedGenotyper -R reference.fasta -I files.list -A AlleleBalance -pnrm EXACT_GENERAL_PLOIDY -ploidy 1 -glm SNP -o files.vcf


java -jar GenomeAnalysisTK.jar -T VariantFiltration -R reference.fasta -V files.vcf --filterExpression "((DP-MQ0)<10) || ((MQ0/(1.0*DP))>=0.8) || (ABHom <0.8) || (Dels >0.5) || (QUAL > 90)" --filterName LowConfidence -o files_filtered.vcf


vcftools --vcf files_filtered.vcf --recode --keep-INFO-all


for file in files_filtered.vcf;do for sample in `/home/bcftools/bcftools-git/bcftools view -h $file | grep "^#CHROM" | cut  -f10-`; do /home/bcftools/bcftools-git/bcftools view -c1 -s $sample -o ${file/.vcf*/.$sample.vcf} $file;done;done

vcf_filter_module.py 9 in.vcf out.vcf

java -jar GenomeAnalysisTK.jar -R reference.fasta -T CombineVariants –V vcf1 vcf2 vcf3 -genotypeMergeOptions UNIQUIFY –o master.vcf

zcat master.vcf.gz | vcf-to-tab > snps.tab

vcf_tab_to_fasta_alignment.pl -i snps.tab > all_snps.fasta
