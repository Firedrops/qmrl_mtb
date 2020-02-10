#!/usr/bin/env nextflow
/*call with nextflow run main.nf -with-docker qimr_mtb:slim */

/**********
 * Define the default parameters
 **********/

params.genome     = "$baseDir/data/H37Rv_refe.fasta"
params.reads      = "$baseDir/data/*_R{1,2}.fastq.gz"
params.results    = "results"
params.gatk       = '/usr/local/bin/GenomeAnalysisTK.jar'
params.gatk_launch = "java -jar $params.gatk"
params.output  = "outdir"
TRIMMOMATIC='/sw/trimmomatic/Trimmomatic-0.36/'
TRIMMOMATIC='/trimmomatic/classes/'
threads=1
adapters="/sw/trimmomatic/Trimmomatic-0.36/adapters/NexteraPE-PE.fa"
adapters="/trimmomatic/adapters/NexteraPE-PE.fa"
leading=10
trailing=10
minlen=36

log.info """\
C A L L I N G S  -  N F    v 1.0
================================
genome   : $params.genome
reads    : $params.reads
variants : $params.variants
blacklist: $params.blacklist
results  : $params.results
gatk     : $params.gatk
"""

/**********
 *  Parse the input parameters
 **********/

GATK            = params.gatk_launch
genome_file     = file(params.genome)

Channel
    .fromFilePairs( params.reads, flat: true )
    .ifEmpty { exit 1, "Read pair files could not be found: ${params.reads}" }
    .set { reads }


/**********
 * PART 1: Data preparation
 *
 * Process 1A: Create a FASTA genome index (.fai) with samtools for GATK
 *  samtools faidx /data/reference.fasta
 **********/

process '1A_prepare_genome_samtools' {
  tag "$genome.baseName"

  input:
      file genome from genome_file

  output:
      file "${genome}.fai" into genome_index_ch

  script:
  """
  samtools faidx ${genome}
  """
}



/**********
 * Process 1B: Create a bwa index (.bwt) with bwa
 *  bwa index /data/reference.fasta
 *
 **********/

process '1B_prepare_genome_bwa' {
  tag "$genome.baseName"

  input:
      file genome from genome_file

  output:
      file "${genome}.bwt" into genome_index_bwt_ch

  script:
  """
  bwa index  ${genome}
  """
}

/**********
 * Process 1C: Create a picard index (.dict) with picard
 *  picard CreateSequenceDictionary R=/data/reference.fasta O=/data/reference.dict
 *
 **********/
process '1C_prepare_genome_picard' {
  tag "$genome.baseName"
  label 'mem_xlarge'

  input:
      file genome from genome_file
  output:
      file "${genome.baseName}.dict" into genome_dict_ch

  script:
  """
  picard CreateSequenceDictionary R=$genome O=${genome.baseName}.dict
  """
}


/*
 *  END OF PART 1
 **********/


/**********
 * PART 2: Alignment
 *
 **********/

/**********
 * SKIPPED Process 2A: trimmomatic SKIPPED
 * from https://github.com/cdeanj/amrplusplus/blob/master/main.nf
 *java -jar /trimmomatic/classes/trimmomatic.jar PE -phred33
 *-trimlog ${file}_log.txt ${file}_R1.fastq.gz ${file}_R2.fastq.gz
 *  ${file}_paired_R1.fastq.gz ${file}_unpaired_R1.fastq.gz ${file}_paired_R2.fastq.gz ${file}_unpaired_R2.fastq.gz
 * ILLUMINACLIP:/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36
 *
 *
 * process RunQC {
 *    tag { sample_id }
 *
 *    publishDir "${params.output}/RunQC", mode: 'copy', pattern: '*.fastq',
 *        saveAs: { filename ->
 *            if(filename.indexOf("P.fastq") > 0) "Paired/$filename"
 *            else if(filename.indexOf("U.fastq") > 0) "Unpaired/$filename"
 *            else {}
 *        }
 *
 *    input:
 *        set sample_id, file(forward), file(reverse) from reads
 *
 *    output:
 *        set sample_id, file("${sample_id}.1P.fastq"), file("${sample_id}.2P.fastq") into (paired_fastq)
 *        set sample_id, file("${sample_id}.1U.fastq"), file("${sample_id}.2U.fastq") into (unpaired_fastq)
 *        file("${sample_id}.trimmomatic.stats.log") into (trimmomatic_stats)
 *
 *    """
 *    java -jar ${TRIMMOMATIC}/trimmomatic.jar \
 *      PE -phred33 \
 *      $forward $reverse  \
 *      ILLUMINACLIP:${adapters}:2:30:10 \
 *      LEADING:${leading} \
 *      TRAILING:${trailing} \
 *      SLIDINGWINDOW:4:15 \
 *      MINLEN:${minlen} \
 *      2> ${sample_id}.trimmomatic.stats.log
 *    mv ${sample_id}_1P ${sample_id}.1P.fastq
 *    mv ${sample_id}_2P ${sample_id}.2P.fastq
 *    mv ${sample_id}_1U ${sample_id}.1U.fastq
 *    mv ${sample_id}_2U ${sample_id}.2U.fastq
 *    """
 *}
 *
 *trimmomatic_stats.toSortedList().set { trim_stats }
 */

/**********
 * Process 2A bwa mem
 *   bwa mem -t 4 -M -R "@RG\tID:${NAME}\tSM:${NAME}\tPL:Illumina\tLB:001\tPU:001"
 *   ${reference}.fasta ${indir}${NAME}_R1.fastq.gz ${indir}${NAME}_R2.fastq.gz >
 *   ${tempdir}${NAME}.sam
 **********/

 params.pair1 = "$baseDir/data/*_R1.fastq.gz"
 params.pair2 = $baseDir/data/*_R2.fastq.gz"

process '2A_bwa_mem_alignment' {
  tag "$genome.baseName"

     cpus params.cpus.mapping

     input:
     file genome from genome_file, params.pair1, params.pair2

     output:
     file "${genome.baseName}" into bam_files

     script:
     """
     bwa mem -M -R '@RG\tID:${genome.baseName}\tSM:${genome.baseName}' -t ${params.cpus.mapping} ${genome_file} ${pair1} ${pair2} | samtools view -bS -F 4 ${genome.baseName}.bam
     """
}

  /**********
  * Process 2B picard SortSam
  *   picard SortSam I=${tempdir}${NAME}.bam O=${tempdir}${NAME}_sorted.bam SO=coordinate
  **********/

  process '2B_picard_SortSam' {
    tag "$genome.baseName"

       input:
       file "${genome.baseName}.bam" from bam_files

       output:
       file "${genome.baseName}_sorted.bam" into bam_files

       script:
       """
       picard SortSam I=${Genome.baseName}.bam O=${genome.baseName}_sorted.bam SO=coordinate
       """
  }

  /**********
  * Process 2C picard MarkDuplicates
  *   picard MarkDuplicates I=${tempdir}${NAME}_sorted.bam O=${tempdir}${NAME}_dup.bam M=${tempdir}${NAME}.txt
  **********/

  process '2C_picard_MarkDuplicates' {
    tag "$genome.baseName"

       input:
       file "${genome.baseName}_sorted.bam" from sortedbam_files

       output:
       file "${genome.baseName}_dup.bam" into dupbam_files, file "${genome.baseName}_metric.txt" into mectric_files

       script:
       """
       picard MarkDuplicates I=${genome.baseName}_sorted.bam O=${genome.baseName}_dup.bam M=${genome.baseName}_metric.txt
       """
  }

  /**********
  * Process 2C picard BuildBamIndex
  *   picard BuildBamIndex I=${tempdir}${NAME}_dup.bam
  **********/

  process '2C_picard_MarkDuplicates' {
    tag "$genome.baseName"

       input:
       file "${genome.baseName}_dup.bam" from dupbam_files

       output:
       file "${genome.baseName}.bai" into bai_files

       script:
       """
       picard BuildBamIndex I=${genome.baseName}_dup.bam O=${genome.baseName}.bai
       """
  }

  /**********
   * PART 3: Filtration
   *
   **********/

  /**********
  * Process 3A gatk RealignerTarget
  *   gatk3 -T RealignerTargetCreator -R ${reference}.fasta -I ${tempdir}${NAME}_dup.bam -o ${tempdir}${NAME}.intervals
  **********/

  process '3A_gatk_RealignerTarget' {
    tag "$genome.baseName"

       input:
       file "${genome.baseName}_dup.bam" from dupbam_files, $genome

       output:
       file "${genome.baseName}.intervals" into interval_files

       script:
       """
       gatk3 -T RealignerTargetCreator -R $genome -I ${genome.baseName}_dup.bam -o ${genome.baseName}.intervals
       """
  }

  /**********
  * Process 3B gatk InDelRealigner
  *   gatk3 -T IndelRealigner -R ${reference}.fasta -I ${tempdir}${NAME}_dup.bam -targetIntervals ${tempdir}${NAME}.intervals -o ${outdir}${NAME}_dup_alig.bam
  **********/

  process '3B_gatk_InDelRealigner' {
    tag "$genome.baseName"

       input:
       file "${genome.baseName}_dup.bam" from dupbam_files, $genome, file "${genome.baseName}.intervals" from interval_files

       output:
       file "${genome.baseName}_dup_alig.bam" into dup_alig_files

       script:
       """
       gatk3 -T IndelRealigner -R $genome -I ${genome.baseName}_dup.bam -targetIntervals ${genome.baseName}.intervals -o ${genome.baseName}_dup_alig.bam
       """
  }


    /**********
    * Process 3C gatk UnifiedGenotyper
    *   gatk3 -T UnifiedGenotyper -R ${reference}.fasta -I ${outdir}${NAME}_dup_alig.bam -A AlleleBalance -pnrm EXACT_GENERAL_PLOIDY -ploidy 1 -glm SNP -o ${outdir}${NAME}.vcf
    **********/
    process '3D_gatk_UnifiedGenotyper' {
      tag "$genome.baseName"

         input:
         file "${genome.baseName}_dup_alig.bam" from dup_alig_files, $genome

         output:
         file "${genome.baseName}.vcf" into vcf_files

         script:
         """
         gatk3 -T UnifiedGenotyper -R $genome -I ${genome.baseName}_dup_alig.bam -A AlleleBalance -pnrm EXACT_GENERAL_PLOIDY -ploidy 1 -glm SNP -o ${outdir}${genome.baseName}.vcf
         """
    }

    /**********
    * Process 3D gatk VariantFiltration
    *   gatk3 -T VariantFiltration -R ${reference}.fasta -V ${outdir}${NAME}.vcf --filterExpression "((DP-MQ0)<10) || ((MQ0/(1.0*DP))>=0.8) || (ABHom <0.8) || (Dels >0.5) || (QUAL > 90)" --filterName LowConfidence -o ${outdir}${NAME}_filtered.vcf
    **********/
    process '3D_gatk_UnifiedGenotyper' {
      tag "$genome.baseName"

         input:
         file "${genome.baseName}.vcf" from vcf_files, $genome

         output:
         file "${genome.baseName}_filtered.vcf" into vcf_filtered_files

         script:
         """
         gatk3 -T VariantFiltration -R $genome -V ${genome.baseName}.vcf --filterExpression "((DP-MQ0)<10) || ((MQ0/(1.0*DP))>=0.8) || (ABHom <0.8) || (Dels >0.5) || (QUAL > 90)" --filterName LowConfidence -o $(genome.baseName}_filtered.vcf
    }

  /**********
  * PART 4: Post-Processing
  *
  **********/

     /**********
     * Process 4A vcftools recode
     *   vcftools --vcf ${outdir}${NAME}_filtered.vcf --recode --keep-INFO-all
     **********/
     process '4A_vcf_recode' {
       tag "$genome.baseName"

          input:
          file "${genome.baseName}_filtered.vcf" from vcf_filtered_files

          output:
          file "${genome.baseName}.recode.vcf" into vcf_recode_files

          script:
          """
          vcftools --vcf ${genome.baseName}_filtered.vcf --out ${genome.baseName}.recode.vcf --recode --keep-INFO-all
          """
     }

     /**********
     * Process 4B bcftools view
     *   bcftools view -c1 -s ${NAME} -o ${tempdir}${NAME}_in.vcf ${tempdir}${NAME}_out.recode.vcf
     **********/
     process '4B_vcf_recode' {
       tag "$genome.baseName"

          input:
          file "${genome.baseName}" from vcf_filtered_files

          output:
          file "${genome.baseName}_in.vcf" into vcf_in_files, file "${tempdir}${NAME}_out.vcf" into vcf_out_files

          script:
          """
          bcftools view -c1 -s {genome.baseName}-o ${genome.baseName}_in.vcf ${tempdir}${NAME}_out.vcf
          """
     }

     /**********
     * Process 4C python filter
     *   python vcf_filter_module.py 9 ${tempdir}${NAME}_in.vcf ${tempdir}${NAME}_out.vcf
     **********/
     process '4C_python_filter' {
       tag "$genome.baseName"

          input:
          file "${genome.baseName}_in.vcf" from vcf_in_files, file "${genome.baseName}_out.vcf" from vcf_out_files

          output:
          I DON'T KNOW

          script:
          """
          python vcf_filter_module.py 9 ${genome.baseName}_in.vcf ${genome.baseName}_out.vcf
          """
     }

     /**********
      * PART 5: Consolidation
      *
      **********/

     /**********
     * Process 5A bgzip
     *   cat ${tempdir}${NAME}_out.vcf | bgzip -c > ${outdir}${NAME}_master.vcf.gz
     **********/

     process '5A_bgzp_vcf' {
       tag "$genome.baseName"

          input:
          file "${genome.baseName}_out.vcf" from vcf_out_files

          output:
          file "${genome.baseName}_master.vcf.gz" into vcf_master_files

          script:
          """
          cat ${genome.baseName}_out.vcf | bgzip -c > ${genome.baseName}_master.vcf.gz
          """
     }

     /**********
     * Process 5B tab conversion
     *   zcat ${outdir}${NAME}_master.vcf.gz | vcf-to-tab > ${tempdir}${NAME}_snps.tab
     **********/

     process '5B_vcf_to_tab' {
       tag "$genome.baseName"

          input:
           file "${genome.baseName}_master.vcf.gz" from vcf_master_files

          output:
          file "${genome.baseName}_snps.tab" into snps_files

          script:
          """
          zcat ${genome.baseName}_master.vcf.gz | vcf-to-tab > ${genome.baseName}_snps.tab
          """
     }
     

     /**********
     * Process 5C all snps fasta
     *   ./vcf_tab_to_fasta_alignment.pl -i ${tempdir}${NAME}_snps.tab > ${outdir}${NAME}_all_snps.fasta
     **********/

     process '5B_vcf_to_tab' {
       tag "$genome.baseName"

          input:
          file "${genome.baseName}_snps.tab" from snps_files

          output:
          file "${genome.baseName}_all_snps.fasta" into all_snps_files

          script:
          """
          ./vcf_tab_to_fasta_alignment.pl -i ${genome.baseName}_snps.tab > ${genome.baseName}_all_snps.fasta
          """
     }



workflow.onComplete { 
	println ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
