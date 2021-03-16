#!/usr/bin/env nextflow

///////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Nextflow Pipeline for HybPiper  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////



def helpMessage() {
    log.info """

    Usage:
    The typical command for running the pipeline is as follows:

    nextflow run hybpiper_pipeline_v1_7_NO_INTRONERATE.nf -c hybpiper.config --illumina_reads_directory <directory> --target_file <file>
    -profile <profile>

    Mandatory arguments:

      ##############################################################################################

      --illumina_reads_directory <directory>    Path to folder containing illumina read file(s)
      --target_file <file>                      File containing fasta sequences of target genes

      ##############################################################################################

    Optional arguments:

      -profile <profile>                        Configuration profile to use. Can use multiple (comma separated)
                                                Available: standard (default), slurm

      --cleanup                                 Run 'cleanup.py' for each gene directory after 'reads_first.py'

      --nosupercontigs                          Do not create supercontigs. Use longest Exonerate hit only. Default is off.  

      --memory <int>                            Memory (RAM) amount in GB to use for bbmap.sh with exonerate_hits.py. Default is 1 GB
      
      --discordant_reads_edit_distance <int>    Minimum number of base differences between one read of a read pair vs the supercontig reference for a read pair to be flagged as discordant. Default is 5
      
      --discordant_reads_cutoff <int>           Minimum number of discordant reads pairs required to flag a supercontigs as a potential hybrid of contigs from multiple paralogs. Default is 5

      --merged                                  Merge forward and reverse reads, and run SPAdes assembly with merged and unmerged (the latter in interleaved format) data. Default is off

      --paired_and_single                       Use when providing both paired R1 and R2 read files as well as a file of single-end reads for each sample       

      --single_only                             Use when providing providing only a folder of single-end reads 

      --outdir <directory_name>                 Specify the name of the pipeline results directory. Default is 'results'                                 

      --read_pairs_pattern <pattern>            Provide a comma-separated read pair pattern for matching fowards and reverse paired-end readfiles. Default is 'R1,R2'

      --single_pattern <pattern>                Provide a pattern for matching single-end read files. Default is 'single'

      --use_blastx                              Use a protein target file and map reads to targets with BLASTx. Default is a nucleotide target file and mapping of reads to targets using BWA

      --num_forks <int>                         Specify the number of parallel processes (e.g. concurrent runs of reads.first.py) to run at any one time. Can be used to prevent Nextflow from using all the threads/cpus on your machine. Default is to use the maximum number possible      

      --cov_cutoff <int>                        Coverage cutoff to pass to the SPAdes assembler. Default is 8

      --blastx_evalue <value>                   Evalue to pass to blastx when using blastx mapping (i.e. when the --use_blastx flag is specified). Default is 1e-4

      --paralog_warning_min_len_percent <decimal> 
                                                Minimum length percentage of a contig vs reference protein length for a paralog warning to be generated and a putative paralog contig to be recovered. Default is 0.75 

      --translate_target_file_for_blastx        Translate a nucleotide target file. If set, the --use_blastx is set by default. Default is off

      --use_trimmomatic                         Trim forwards and reverse reads using Trimmomatic. Default is off

      --trimmomatic_leading_quality <int>       Cut bases off the start of a read, if below this threshold quality. Default is 3

      --trimmomatic_trailing_quality <int>      Cut bases off the end of a read, if below this threshold quality. Default 3

      --trimmomatic_min_length <int>            Drop a read if it is below this specified length. Default is 36

      --trimmomatic_sliding_window_size <int>   Size of the sliding window used by Trimmomatic; specifies the number of bases to average across. Default is 4

      --trimmomatic_sliding_window_quality <int>
                                                Specifies the average quality required within the sliding window. Default is 20

      --run_intronerate                        Run intronerate.py to recover (hopefully) intron and supercontig  sequences. Default is off, and so results `subfolders 09_sequences_intron` and `10_sequences_supercontig` will be empty


    """.stripIndent()
}


/**
* @function printAllMethods
* @purpose Prints an objects class name and then list the associated class functions.
* From https://bateru.com/news/2011/11/code-of-the-day-groovy-print-all-methods-of-an-object/
**/
// Filename: printAllMethodsExample.groovy
void printAllMethods( obj ){
    if( !obj ){
    println( "Object is null\r\n" );
    return;
    }
  if( !obj.metaClass && obj.getClass() ){
        printAllMethods( obj.getClass() );
    return;
    }
  def str = "class ${obj.getClass().name} functions:\r\n";
  obj.metaClass.methods.name.unique().each{ 
    str += it+"(); "; 
  }
  println "${str}\r\n";
}


params.help = false
params.illumina_reads_directory = false
params.target_file = false

// Check that input directories are provided
if (params.help || !params.illumina_reads_directory || !params.target_file) {
  helpMessage()
  exit 0
}

// Check that paralog_warning_min_len_percent value is a decimal between 0 and 1
if (params.paralog_warning_min_len_percent < 0 || params.paralog_warning_min_len_percent >1) {
println("The value for --paralog_warning_min_len_percent should be between 0 and 1. Your value is ${params.paralog_warning_min_len_percent}")
exit 0
}

// Check that non-overlapping options are provided
if (params.single_only && params.paired_and_single) {
  println('Please use --single_only OR --paired_and_single, not both!')
  exit 0
}

// Don't allow params.paired_and_single and params.use_trimmomatic
if (params.paired_and_single && params.use_trimmomatic) {
  println("Trimmomatic can't be used with paired plus single reads yet - let me know if this would be useful!")
  exit 0
}

// Check for unrecognised pararmeters
allowed_params = ["cleanup", "nosupercontigs", "memory","discordant_reads_edit_distance", "discordant_reads_cutoff", "merged", "paired_and_single", "single_only", "outdir", "illumina_reads_directory", "target_file", "help", "memory", "read_pairs_pattern", "single_pattern", "use_blastx", "num_forks", "cov_cutoff", "blastx_evalue", "paralog_warning_min_len_percent", "translate_target_file_for_blastx", "use_trimmomatic", "trimmomatic_leading_quality", "trimmomatic_trailing_quality", "trimmomatic_min_length", "trimmomatic_sliding_window_size", "trimmomatic_sliding_window_quality", "run_intronerate"]

params.each { entry ->
  if (! allowed_params.contains(entry.key)) {
      println("The parameter <${entry.key}> is not known");
      exit 0;
  }
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Set up channels for input data  ///////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Target gene sequences file
///////////////////////////////////////////////////////////////////////////////////////////////////////////

Channel
  .fromPath("${params.target_file}", checkIfExists: true)
  .first()
  .set { target_file_ch }

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Corrected Illumina reads
///////////////////////////////////////////////////////////////////////////////////////////////////////////

// Function for grouping reads from multiple lanes, based on shared filename prefix preceeding the first 
// underscore
def getLibraryId( prefix ){
  prefix.split("_")[0]
}

// Unpaired reads only
if (params.single_only) {
  Channel
  .fromPath("${params.illumina_reads_directory}/*_{$params.single_pattern}*.{fastq.gz,fastq,fq.gz,fq}", \
    checkIfExists: true)
  .map { file -> tuple(file.baseName.split('_')[0], file) }
  .groupTuple(sort:true)
  // .view()
  .set { illumina_reads_single_only_ch }
} else {
  illumina_reads_single_only_ch = Channel.empty()
}


// Paired reads and a file of unpaired reads
if (params.paired_and_single) {
  Channel
  .fromFilePairs("${params.illumina_reads_directory}/*_{$params.read_pairs_pattern,$params.single_pattern}*.{fastq.gz,fastq,fq.gz,fq}", flat : true,
  checkIfExists: true, size: 3)
  .set { illumina_paired_reads_with_unpaired_ch }
} else {
  illumina_paired_reads_with_unpaired_ch = Channel.empty()
}


// Paired read only
// Gather the pairs of R1/R2 according to sample ID
if (!params.paired_and_single && !params.single_only) {
Channel
     .fromFilePairs("${params.illumina_reads_directory}/*_{$params.read_pairs_pattern}*.{fastq.gz,fastq,fq.gz,fq}", \
    flat : true, checkIfExists: true)
     .map { prefix, file1, file2 -> tuple(getLibraryId(prefix), file1, file2) }
     .groupTuple(sort:true)
     .set{ illumina_paired_reads_ch_1 }
} else {
  illumina_paired_reads_ch_1 = Channel.empty()
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Create 'namelist.txt' file and associated channel
///////////////////////////////////////////////////////////////////////////////////////////////////////////

if (!params.single_only) {
  Channel
  .fromFilePairs("${params.illumina_reads_directory}/*_{$params.read_pairs_pattern}*.{fastq.gz,fastq,fq.gz,fq}", \
    flat : true, checkIfExists: true)
  .map { prefix, file1, file2 -> tuple(getLibraryId(prefix), file1, file2) }
  .groupTuple(sort:true)
  .collectFile(name: "${params.outdir}/01_namelist/namelist.txt") { item -> item[0] + "\n" }
  .first()
  .set { namelist_ch }
} else {
  Channel
  .fromPath("${params.illumina_reads_directory}/*_{$params.single_pattern}*.{fastq.gz,fastq,fq.gz,fq}", \
    checkIfExists: true)
  .map { file -> tuple(file.baseName.split('_')[0], file) }
  .groupTuple(sort:true)
  .collectFile(name: "${params.outdir}/01_namelist/namelist.txt") { item -> item[0] + "\n" }
  .first()
  .set { namelist_ch }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Channel of gene names for 'paralog_retriever.py' script
///////////////////////////////////////////////////////////////////////////////////////////////////////////

Channel
.fromPath("${params.target_file}", checkIfExists: true)
.splitFasta( record: [id: true, seqString: true ])
.map { it.id.replaceFirst(~/.*-/, '') }
.unique()
.set { gene_names_ch }
// gene_names_ch.view { "value: $it" }


///////////////////////////////////////////////////////////////////////////////////////////////////////////
// If the flag '--translate_target_file_for_blastx' is set, translate nucleotide target file  
///////////////////////////////////////////////////////////////////////////////////////////////////////////

process translate_target_file {
 echo true
 label 'in_container'
 publishDir "${params.outdir}/00_translated_target_file", mode: 'copy'

 when:
 params.translate_target_file_for_blastx

 input:
 file(target_file_nucleotides) from target_file_ch

 output:
 file("target_file_translated.fasta") into target_file_translated_ch
 file("translation_warnings.txt")

 script:
 """
 #!/usr/bin/env python

from Bio import SeqIO
 
translated_seqs_to_write = []
with open("${target_file_nucleotides}", 'r') as target_file_nucleotides:
  seqs = SeqIO.parse(target_file_nucleotides, 'fasta')
  with open('translation_warnings.txt', 'w') as translation_warnings:
    for seq in seqs:
      if len(seq.seq) % 3 != 0:
        translation_warnings.write(f"WARNING: sequence for gene {seq.name} is not a multiple of 3. Translating anyway...\\n")
      protein_translation = seq.translate()
      protein_translation.name = seq.name
      protein_translation.id = seq.id
      protein_translation.description = 'translated sequence from nucleotide target file'
      num_stop_codons = protein_translation.seq.count('*')
      if num_stop_codons != 0:
        translation_warnings.write(f'WARNING: stop codons present in translation of sequence {seq.name}, please check\\n')
      translated_seqs_to_write.append(protein_translation)

with open('target_file_translated.fasta', 'w') as translated_handle:
  SeqIO.write(translated_seqs_to_write, translated_handle, 'fasta')

"""
}

// set correct target_file_ch based on whether using BWA or BLASTx
if (!params.translate_target_file_for_blastx) {
  target_file_ch = target_file_ch
} else {
  target_file_ch = target_file_translated_ch
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Concatenate sample files from multiple lanes if necessary
///////////////////////////////////////////////////////////////////////////////////////////////////////////

//Paired R1 and R2 reads only:
process combine_lanes_paired_only {
  label 'in_container'
  // echo true
  publishDir "$params.outdir/02_reads_combined_lanes", mode: 'copy', pattern: "*.fastq*"

  if (params.num_forks) {
      maxForks params.num_forks
  }

  input:
  tuple val(prefix), file(reads_R1), file(reads_R2) from illumina_paired_reads_ch_1

  output:
  tuple val(prefix), file("*R1.fastq*"), file("*R2.fastq*") into combined_lane_paired_reads_ch_1, combined_lane_paired_reads_ch_2

  script:
  """
  first_file=\$(echo $reads_R1 | cut -d' ' -f1)
  echo \$first_file

  if [[ \$first_file = *.gz ]]
    then 
      cat $reads_R1 > ${prefix}_combinedLanes_R1.fastq.gz
      cat $reads_R2 > ${prefix}_combinedLanes_R2.fastq.gz
  fi

  if [[ \$first_file = *.fq ]] || [[ \$first_file = *.fastq ]]
    then 
      cat $reads_R1 > ${prefix}_combinedLanes_R1.fastq
      cat $reads_R2 > ${prefix}_combinedLanes_R2.fastq
  fi
  """
}


//Single reads only:
process combine_lanes_single_only {
  label 'in_container'
  // echo true
  publishDir "$params.outdir/02_reads_combined_lanes", mode: 'copy', pattern: "*.fastq*"

  if (params.num_forks) {
      maxForks params.num_forks
  }

  input:
  tuple val(prefix), file(reads_single) from illumina_reads_single_only_ch

  output:
  tuple val(prefix), file("*single.fastq*")  into combined_lane_single_reads_ch_1, combined_lane_single_reads_ch_2

  script:
  """
  first_file=\$(echo $reads_single | cut -d' ' -f1)
  echo \$first_file

  if [[ \$first_file = *.gz ]]
    then 
      cat $reads_single > ${prefix}_combinedLanes_single.fastq.gz
  fi

  if [[ \$first_file = *.fq ]] || [[ \$first_file = *.fastq ]]
    then 
      cat $reads_single > ${prefix}_combinedLanes_single.fastq
  fi
  """
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Run optional Trimmomatic step
///////////////////////////////////////////////////////////////////////////////////////////////////////////

// R1 and R2 reads
process trimmomatic_paired {
  // echo true
  label 'in_container'
  publishDir "$params.outdir/03b_trimmomatic_paired_and_single_reads", mode: 'copy', pattern: "*_paired.fq*"
  publishDir "$params.outdir/03b_trimmomatic_paired_and_single_reads", mode: 'copy', pattern: "*_both_unpaired.fq*"
  publishDir "$params.outdir/03c_trimmomatic_logs", mode: 'copy', pattern: "*.log"

  if (params.num_forks) {
      maxForks params.num_forks
  }

  when:
  params.use_trimmomatic

  input:
  tuple val(prefix), file(combined_reads_R1), file(combined_reads_R2) from combined_lane_paired_reads_ch_1

  output:
  file("*")
  tuple val(prefix), file("*forward_paired*"), file("*reverse_paired*"), file("*both_unpaired*") into trimmed_paired_and_orphaned_ch

  script:
  """
  R1=${combined_reads_R1}
  R2=${combined_reads_R2}

  if [[ \$R1 = *.gz ]]
    then 
      sampleID=\${R1%_R[1,2].fastq.gz}
      output_forward_paired=\${R1%.fastq.gz}_forward_paired.fq.gz
      output_reverse_paired=\${R2%.fastq.gz}_reverse_paired.fq.gz
      output_forward_unpaired=\${R1%.fastq.gz}_forward_unpaired.fq.gz
      output_reverse_unpaired=\${R2%.fastq.gz}_reverse_unpaired.fq.gz
      output_both_unpaired=\${sampleID}_both_unpaired.fq.gz
    else
      sampleID=\${R1%_R[1,2].fastq}
      output_forward_paired=\${R1%.fastq}_forward_paired.fq
      output_reverse_paired=\${R2%.fastq}_reverse_paired.fq
      output_forward_unpaired=\${R1%.fastq}_forward_unpaired.fq
      output_reverse_unpaired=\${R2%.fastq}_reverse_unpaired.fq
      output_both_unpaired=\${sampleID}_both_unpaired.fq
  fi

  echo -e ">PrefixPE/1\nTACACTCTTTCCCTACACGACGCTCTTCCGATCT\n>PrefixPE/2\nGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT\n\
  >PE1\nTACACTCTTTCCCTACACGACGCTCTTCCGATCT\n>PE1_rc\nAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA\n\
  >PE2\nGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT\n>PE2_rc\nAGATCGGAAGAGCACACGTCTGAACTCCAGTCA" > TruSeq3-PE-2.fa
  trimmomatic PE -phred33 -threads ${task.cpus} \
  ${combined_reads_R1} ${combined_reads_R2} \${output_forward_paired} \${output_forward_unpaired} \
  \${output_reverse_paired} \${output_reverse_unpaired} \
  ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:1:true LEADING:${params.trimmomatic_leading_quality} TRAILING:${params.trimmomatic_trailing_quality} SLIDINGWINDOW:${params.trimmomatic_sliding_window_size}:${params.trimmomatic_sliding_window_quality} MINLEN:${params.trimmomatic_min_length} 2>&1 | tee \${sampleID}.log 
  cat \${output_forward_unpaired} \${output_reverse_unpaired} > \${output_both_unpaired}
  """
}


// Single reads only
process trimmomatic_single {
  // echo true
  label 'in_container'
  publishDir "$params.outdir/03a_trimmomatic_single_reads", mode: 'copy', pattern: "*_single*"
  publishDir "$params.outdir/03b_trimmomatic_logs", mode: 'copy', pattern: "*.log"

  if (params.num_forks) {
      maxForks params.num_forks
  }

  when:
  params.use_trimmomatic

  input:
  tuple val(prefix), file(combined_reads_single) from combined_lane_single_reads_ch_1

  output:
  file("*")
  tuple val(prefix), file("*single*") into trimmed_single_ch

  script:
  """
  single=${combined_reads_single}


  if [[ \$single = *.gz ]]
    then 
      output_single=${prefix}_trimmed_single.fastq.gz
    else
      output_single=${prefix}_trimmed_single.fastq
  fi

  echo -e ">TruSeq3_IndexedAdapter\nAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC\n>TruSeq3_UniversalAdapter\nAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA\n" > TruSeq3-SE.fa
  trimmomatic SE -phred33 -threads ${task.cpus} \
  ${combined_reads_single} \${output_single} ILLUMINACLIP:TruSeq3-SE.fa:2:30:10:1:true LEADING:${params.trimmomatic_leading_quality} TRAILING:${params.trimmomatic_trailing_quality} SLIDINGWINDOW:${params.trimmomatic_sliding_window_size}:${params.trimmomatic_sliding_window_quality} MINLEN:${params.trimmomatic_min_length} 2>&1 | tee ${prefix}.log 
  """
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Set up conditional channels for single end reads (trimmed vs non-trimmed)
///////////////////////////////////////////////////////////////////////////////////////////////////////////

(single_channel_1, single_channel_2) = (params.use_trimmomatic ? 
  [Channel.empty(), trimmed_single_ch] : [combined_lane_single_reads_ch_2, Channel.empty()] )



///////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////  Run HybPiper  ///////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// reads_first.py
///////////////////////////////////////////////////////////////////////////////////////////////////////////

// Single-end reads only:
process reads_first_with_single_end_only {
  // echo true
  label 'in_container'
  publishDir "${params.outdir}/06_summary_stats", mode: 'copy', pattern: "${pair_id}/${pair_id}_genes_with_supercontigs.csv"
  publishDir "${params.outdir}/06_summary_stats", mode: 'copy', pattern: "${pair_id}/${pair_id}_supercontigs_with_discordant_reads.csv"

  if (params.num_forks) {
      maxForks params.num_forks
  }

  when:
  params.single_only

  input:
  file(target_file) from target_file_ch
  tuple val(prefix), file(reads_single) from single_channel_1.mix(single_channel_2)

  output:
  file("${prefix}") optional true into (reads_first_with_single_only_ch_1, reads_first_with_single_only_ch_2, reads_first_with_single_only_ch_3, reads_first_with_single_only_ch_4)
  file("${prefix}/${pair_id}_genes_with_supercontigs.csv") optional true
  file("${prefix}/${pair_id}_supercontigs_with_discordant_reads.csv") optional true

  script:

  def command_list = []
  if (params.nosupercontigs) {
     command_list << "--nosupercontigs"
     }
  if (params.memory) {
     command_list << "--memory ${params.memory}"
     }
  if (params.discordant_reads_edit_distance) {
     command_list << "--discordant_reads_edit_distance ${params.discordant_reads_edit_distance}"
     }
  if (params.discordant_reads_cutoff) {
     command_list << "--discordant_reads_cutoff ${params.discordant_reads_cutoff}"
     } 
  if (params.merged) {
     command_list << "--merged"
     }
  if (!params.use_blastx && !params.translate_target_file_for_blastx) {
    command_list << "--bwa"
  }
  if (params.blastx_evalue) {
    command_list << "--evalue ${params.blastx_evalue}"
  }
  if (params.paralog_warning_min_len_percent) {
    command_list << "--paralog_warning_min_length_percentage ${params.paralog_warning_min_len_percent}"
  }
  if (params.cov_cutoff) {
    command_list << "--cov_cutoff ${params.cov_cutoff}"
  }
  if (params.cleanup) {
     cleanup = "cleanup.py ${prefix}"
  } else {
     cleanup = ''
  }
  reads_first_command = "python /HybPiper/reads_first.py -b ${target_file} -r ${reads_single} --prefix ${prefix} --cpu ${task.cpus} " + command_list.join(' ')

  """
  echo ${reads_first_command}
  ${reads_first_command}
  ${cleanup}
  """
 }


// R1 and R2 reads and a file of single-end reads:
process reads_first_with_unpaired {
  //echo true
  label 'in_container'
  publishDir "${params.outdir}/06_summary_stats", mode: 'copy', pattern: "${pair_id}/${pair_id}_genes_with_supercontigs.csv"
  publishDir "${params.outdir}/06_summary_stats", mode: 'copy', pattern: "${pair_id}/${pair_id}_supercontigs_with_discordant_reads.csv"

  if (params.num_forks) {
      maxForks params.num_forks
  }

  when:
  (params.use_trimmomatic || params.paired_and_single)

  input:
  file(target_file) from target_file_ch
  tuple pair_id, file(reads_R1), file(reads_R2), file(reads_unpaired) from illumina_paired_reads_with_unpaired_ch.mix(trimmed_paired_and_orphaned_ch)

  output:
  file("${pair_id}") optional true into (reads_first_with_unPaired_ch_1, reads_first_with_unPaired_ch_2, reads_first_with_unPaired_ch_3, reads_first_with_unPaired_ch_4)
  file("${pair_id}/${pair_id}_genes_with_supercontigs.csv") optional true
  file("${pair_id}/${pair_id}_supercontigs_with_discordant_reads.csv") optional true

  script:
  def command_list = []
  if (params.nosupercontigs) {
     command_list << "--nosupercontigs"
     }
  if (params.memory) {
     command_list << "--memory ${params.memory}"
     }
  if (params.discordant_reads_edit_distance) {
     command_list << "--discordant_reads_edit_distance ${params.discordant_reads_edit_distance}"
     }
  if (params.discordant_reads_cutoff) {
     command_list << "--discordant_reads_cutoff ${params.discordant_reads_cutoff}"
     } 
  if (params.merged) {
     command_list << "--merged"
     }
  if (!params.use_blastx && !params.translate_target_file_for_blastx) {
    command_list << "--bwa"
  }
  if (params.blastx_evalue) {
    command_list << "--evalue ${params.blastx_evalue}"
  }
  if (params.paralog_warning_min_len_percent) {
    command_list << "--paralog_warning_min_length_percentage ${params.paralog_warning_min_len_percent}"
  }
  if (params.cov_cutoff) {
    command_list << "--cov_cutoff ${params.cov_cutoff}"
  }
  if (params.cleanup) {
     cleanup = "cleanup.py ${pair_id}"
  } else {
     cleanup = ''
  }
  reads_first_command = "python /HybPiper/reads_first.py -b ${target_file} -r ${reads_R1} ${reads_R2} --unpaired ${reads_unpaired} --prefix ${pair_id} --cpu ${task.cpus} " + command_list.join(' ')

  script:
  """
  echo ${reads_first_command}
  ${reads_first_command}
  ${cleanup}
  """
 } 


///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Just R1 and R2 reads:
///////////////////////////////////////////////////////////////////////////////////////////////////////////

process reads_first_no_unpaired {
  // echo true
  label 'in_container'
  publishDir "${params.outdir}/06_summary_stats", mode: 'copy', pattern: "${pair_id}/${pair_id}_genes_with_supercontigs.csv"
  publishDir "${params.outdir}/06_summary_stats", mode: 'copy', pattern: "${pair_id}/${pair_id}_supercontigs_with_discordant_reads.csv"

  if (params.num_forks) {
      maxForks params.num_forks
  }

  when:
  (!params.paired_and_single && !params.single_only && !params.use_trimmomatic)

  input:
  file(target_file) from target_file_ch
  tuple pair_id, file(reads_R1), file(reads_R2) from combined_lane_paired_reads_ch_2

  output:
  file("${pair_id}") optional true into (reads_first_ch_1, reads_first_ch_2, reads_first_ch_3)
  file("${pair_id}/${pair_id}_genes_with_supercontigs.csv") optional true
  file("${pair_id}/${pair_id}_supercontigs_with_discordant_reads.csv") optional true

  script:

  def command_list = []
  if (params.nosupercontigs) {
     command_list << "--nosupercontigs"
     }
  if (params.memory) {
     command_list << "--memory ${params.memory}"
     }
  if (params.discordant_reads_edit_distance) {
     command_list << "--discordant_reads_edit_distance ${params.discordant_reads_edit_distance}"
     }
  if (params.discordant_reads_cutoff) {
     command_list << "--discordant_reads_cutoff ${params.discordant_reads_cutoff}"
     }
  if (params.merged) {
     command_list << "--merged"
     }
  if (!params.use_blastx && !params.translate_target_file_for_blastx) {
    command_list << "--bwa"
  }
  if (params.blastx_evalue) {
    command_list << "--evalue ${params.blastx_evalue}"
  }
  if (params.paralog_warning_min_len_percent) {
    command_list << "--paralog_warning_min_length_percentage ${params.paralog_warning_min_len_percent}"
  }
  if (params.cov_cutoff) {
    command_list << "--cov_cutoff ${params.cov_cutoff}"
  }
  if (params.cleanup) {
     cleanup = "cleanup.py ${pair_id}"
  } else {
     cleanup = ''
  }
  reads_first_command = "python /HybPiper/reads_first.py -b ${target_file} -r ${reads_R1} ${reads_R2} --prefix ${pair_id} --cpu ${task.cpus} " + command_list.join(' ')

  script:
  """
  echo "about to try command: ${reads_first_command}"
  ${reads_first_command}
  ${cleanup}
  """
 }


///////////////////////////////////////////////////////////////////////////////////////////////////////////
// get_seq_lengths.py and gene_recovery_heatmap_ggplot.R 
///////////////////////////////////////////////////////////////////////////////////////////////////////////

process visualise {
 // echo true
 label 'in_container'
 publishDir "${params.outdir}/05_visualise", mode: 'copy'

 input:
 file(reads_first) from reads_first_ch_1.collect().mix(reads_first_with_unPaired_ch_1).collect().mix(reads_first_with_single_only_ch_1).collect()
 file(target_file) from target_file_ch
 file(namelist) from namelist_ch


 output:
 file("seq_lengths.txt") into seq_lengths_ch
 file("heatmap.png")

 script:
 """
 python /HybPiper/get_seq_lengths.py ${target_file} ${namelist} dna > seq_lengths.txt
 Rscript /HybPiper/gene_recovery_heatmap_ggplot.R
 """
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
// hybpiper_stats.py
///////////////////////////////////////////////////////////////////////////////////////////////////////////

process summary_stats {
 // echo true
 label 'in_container'
 publishDir "${params.outdir}/06_summary_stats", mode: 'copy'

 input:
 file(reads_first) from reads_first_ch_2.collect().mix(reads_first_with_unPaired_ch_2).collect().mix(reads_first_with_single_only_ch_2).collect()
 file(seq_lengths) from seq_lengths_ch
 file(namelist) from namelist_ch

 output:
 file("stats.txt")

 script:
 if (params.translate_target_file_for_blastx || params.use_blastx) {
  """
  python /HybPiper/hybpiper_stats.py ${seq_lengths} ${namelist} --blastx_adjustment > stats.txt
  """
  } else {
  """
  python /HybPiper/hybpiper_stats.py ${seq_lengths} ${namelist} > stats.txt
  """
 } 
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Set up conditional channels to skip or include intronerate.py
///////////////////////////////////////////////////////////////////////////////////////////////////////////

(reads_first_channel_1, reads_first_channel_2) = (params.run_intronerate ? 
  [Channel.empty(), reads_first_ch_3.mix(reads_first_with_unPaired_ch_3).mix(reads_first_with_single_only_ch_3)] : [reads_first_ch_3.mix(reads_first_with_unPaired_ch_3).mix(reads_first_with_single_only_ch_3), Channel.empty()] )

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// intronerate.py
///////////////////////////////////////////////////////////////////////////////////////////////////////////

process intronerate {
 // echo true
 label 'in_container'

 input:
 file(reads_first) from reads_first_channel_2

 output:
 file(reads_first) optional true into (intronate_ch)

 script:
 """
 echo ${reads_first}
 python /HybPiper/intronerate.py --prefix ${reads_first}
 """
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
// paralog_investigator.py
///////////////////////////////////////////////////////////////////////////////////////////////////////////

process paralogs {
  //echo true
  label 'in_container'
  publishDir "${params.outdir}/04_processed_gene_directories", mode: 'copy'


  if (params.num_forks) {
      maxForks params.num_forks
   }

  input:
  file(intronerate_complete) from intronate_ch.mix(reads_first_channel_1)

  output:
  file(intronerate_complete) optional true into (paralogs_ch_1, paralogs_ch_2)

  script:
  """
  python /HybPiper/paralog_investigator.py ${intronerate_complete}
  """
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
// retrieve_sequences.py
///////////////////////////////////////////////////////////////////////////////////////////////////////////

process retrieve_sequences {
  // echo true
  label 'in_container'
  publishDir "${params.outdir}/07_sequences_dna", mode: 'copy', pattern: "*.FNA"
  publishDir "${params.outdir}/08_sequences_aa", mode: 'copy', pattern: "*.FAA"
  publishDir "${params.outdir}/09_sequences_intron", mode: 'copy', pattern: "*introns.fasta"
  publishDir "${params.outdir}/10_sequences_supercontig", mode: 'copy', pattern: "*supercontig.fasta"

  input:
  file(paralog_complete) from paralogs_ch_1.collect()
  file(target_file) from target_file_ch


  output:
  file("*.FNA")
  file("*.FAA")
  file("*.fasta")

  script:
  """
  python /HybPiper/retrieve_sequences.py ${target_file} . dna
  python /HybPiper/retrieve_sequences.py ${target_file} . aa
  python /HybPiper/retrieve_sequences.py ${target_file} . intron
  python /HybPiper/retrieve_sequences.py ${target_file} . supercontig
  """
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
// paralog_retriever.py
///////////////////////////////////////////////////////////////////////////////////////////////////////////

process paralog_retriever {
  //echo true
  label 'in_container'
  publishDir "${params.outdir}/11_paralogs", mode: 'copy', pattern: "*.paralogs.fasta"
  publishDir "${params.outdir}/12_paralogs_noChimeras", mode: 'copy', pattern: "*.paralogs_noChimeras.fasta"
  publishDir "${params.outdir}/12_paralogs_noChimeras/logs", mode: 'copy', pattern: "*mylog*"

  input:
  file(paralog_complete_list) from paralogs_ch_2.collect()
  file(namelist) from namelist_ch
  val(gene_list) from gene_names_ch.collect()

  output:
  file("*.fasta")
  file("*.mylog*")  

  script:
  assert (gene_list in List)
  list_of_names = gene_list.join(' ') // Note that this is necessary so that the list isn't of the form [4471, 4527, etc]
  """
  for gene_name in ${list_of_names}
  do
    python /HybPiper/paralog_retriever.py ${namelist} \${gene_name} > \${gene_name}.paralogs_noChimeras.fasta 2> \${gene_name}.paralogs.fasta
  done
  """
}


////////////////////////////////////  End of script ///////////////////////////////////////////////////////
