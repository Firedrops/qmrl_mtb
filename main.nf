#!/usr/bin/env nextflow
/*call with nextflow run main.nf -with-docker qimr_slim */



/*
 * Define the default parameters
 */ 

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

/*
 *  Parse the input parameters
 */

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
 */

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
 */

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
 */
process '1C_prepare_genome_picard' {
  tag "$genome.baseName"
  label 'mem_xlarge'

  input:
      file genome from genome_file
  output:
      file "${genome.baseName}.dict" into genome_dict_ch

  script:
  """
  picard CreateSequenceDictionary R= $genome O= ${genome.baseName}.dict
  """
}


/*
 *  END OF PART 1
 *********/


/**********
 * PART 2: Alignment
 *
 * Process 2A: trimmomatic
  * from https://github.com/cdeanj/amrplusplus/blob/master/main.nf
 *java -jar /trimmomatic/classes/trimmomatic.jar PE -phred33 
 *-trimlog ${file}_log.txt ${file}_R1.fastq.gz ${file}_R2.fastq.gz
*  ${file}_paired_R1.fastq.gz ${file}_unpaired_R1.fastq.gz ${file}_paired_R2.fastq.gz ${file}_unpaired_R2.fastq.gz
 * ILLUMINACLIP:/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36
*
 */
 process RunQC {
    tag { sample_id }

    publishDir "${params.output}/RunQC", mode: 'copy', pattern: '*.fastq',
        saveAs: { filename ->
            if(filename.indexOf("P.fastq") > 0) "Paired/$filename"
            else if(filename.indexOf("U.fastq") > 0) "Unpaired/$filename"
            else {}
        }
	
    input:
        set sample_id, file(forward), file(reverse) from reads

    output:
        set sample_id, file("${sample_id}.1P.fastq"), file("${sample_id}.2P.fastq") into (paired_fastq)
        set sample_id, file("${sample_id}.1U.fastq"), file("${sample_id}.2U.fastq") into (unpaired_fastq)
        file("${sample_id}.trimmomatic.stats.log") into (trimmomatic_stats)

    """
    java -jar ${TRIMMOMATIC}/trimmomatic.jar \
      PE -phred33 \
      $forward $reverse  \
      ILLUMINACLIP:${adapters}:2:30:10 \
      LEADING:${leading} \
      TRAILING:${trailing} \
      SLIDINGWINDOW:4:15 \
      MINLEN:${minlen} \
      2> ${sample_id}.trimmomatic.stats.log
    mv ${sample_id}_1P ${sample_id}.1P.fastq
    mv ${sample_id}_2P ${sample_id}.2P.fastq
    mv ${sample_id}_1U ${sample_id}.1U.fastq
    mv ${sample_id}_2U ${sample_id}.2U.fastq
    """
}

trimmomatic_stats.toSortedList().set { trim_stats }
