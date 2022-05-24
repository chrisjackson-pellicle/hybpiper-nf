#!/usr/bin/env nextflow

//////////////////////////////////////
//  Nextflow Pipeline for HybPiper2  // 
//////////////////////////////////////

nextflow.enable.dsl=2

def helpMessage() {
    log.info """

    Usage:
    The typical command for running the pipeline is as follows:

    nextflow run hybpiper-rbgv-pipeline.nf \
    -c hybpiper-rbgv.config \
    --illumina_reads_directory <directory> \
    --targetfile_dna <fasta_file> \
    -profile <profile>

    Mandatory arguments:

      ############################################################################

      --illumina_reads_directory <directory>    
                                  Path to folder containing illumina read file(s)

      AND

      --targetfile_dna <file>    File containing fasta sequences of target genes
                                  (nucleotides)
                                        
      OR

      --targetfile_aa <file>     File containing fasta sequences of target genes
                                  (amino-acids)

      #############################################################################

    Optional arguments:

      -profile <profile>          Configuration profile to use. Can use multiple 
                                  (comma separated). Available: standard (default), 
                                  slurm

      --namelist                  A text file containing sample names. Only these 
                                  samples will be processed, By default, all samples 
                                  in the provided <Illumina_reads_directory> 
                                  directory are processed

      --combine_read_files        Group and concatenate read-files via a common prefix. 
                                  Useful if samples have been run across multiple lanes. 
                                  Default prefix is all text preceding the first 
                                  underscore (_) in read filenames

      --combine_read_files_num_fields <int>     
                                  Number of fields (delimited by an underscore) to use 
                                  for combining read files when using the 
                                  `--combine_read_files` flag. Default is 1

      --num_forks <int>           Specify the number of parallel processes (e.g. 
                                  concurrent runs of 'hybpiper assemble') to run at any 
                                  one time. Can be used to prevent Nextflow from using 
                                  all the threads/cpus on your machine. Default is 
                                  to use the maximum number possible  
      
      --outdir <directory_name>                 
                                  Specify the name of the pipeline results directory. 
                                  Default is 'results'        

      --paired_and_single         Use when providing both paired-end R1 and R2 read 
                                  files as well as a file of single-end reads for each 
                                  sample       

      --single_end                Use when providing providing only a folder of 
                                  single-end reads                         

      --read_pairs_pattern <pattern>            
                                  Provide a comma-separated read-pair pattern for 
                                  matching fowards and reverse paired-end readfiles, 
                                  e.g. '1P,2P'. Default is 'R1,R2'

      --single_pattern <pattern>                
                                  Provide a pattern for matching single-end read 
                                  files. Default is 'single'

      #######################################################################################
      ################################ Trimmomatic options: #################################
      #######################################################################################

      --use_trimmomatic           Trim forwards and reverse reads using Trimmomatic.
                                  Default is off

      --trimmomatic_leading_quality <int>       
                                  Cut bases off the start of a read, if below this 
                                  threshold quality.Default is 3

      --trimmomatic_trailing_quality <int>      
                                  Cut bases off the end of a read, if below this 
                                  threshold quality. Default is 3

      --trimmomatic_min_length <int>            
                                  Drop a read if it is below this specified length. 
                                  Default is 36

      --trimmomatic_sliding_window_size <int>   
                                  Size of the sliding window used by Trimmomatic; 
                                  specifies the number of bases to average across. 
                                  Default is 4

      --trimmomatic_sliding_window_quality <int>
                                  Specifies the average quality required within the 
                                  sliding window. Default is 20


      #######################################################################################
      ############################# hybpiper assemble options: ##############################
      #######################################################################################

      --bwa                       Use BWA to search reads for hits to target. Requires
                                  BWA and a target file that is nucleotides!

      --diamond                   Use DIAMOND instead of BLASTx

      --diamond_sensitivity       Use the provided sensitivity for DIAMOND searches. 
                                  Option are: 'mid-sensitive', 'sensitive', 
                                  'more-sensitive', 'very-sensitive', 'ultra-sensitive'

      --distribute_hi_mem         Distributing and writing reads to individual gene 
                                  directories  will be 40-50 percent faster, but can use 
                                  more memory/RAM with large input files

      --evalue                    Evalue to pass to blastx when using blastx mapping, 
                                  i.e., when the --use_blastx or 
                                  --translate_target_file_for_blastx flag is specified. 
                                  Default is 1e-4

      --max_target_seqs           Max target seqs to save in BLASTx search, default is 10

      --cov_cutoff <int>          Coverage cutoff to pass to the SPAdes assembler. 
                                  Default is 8

      --single_cell_assembly      Run SPAdes assemblies in MDA (single-cell) mode. 
                                  Default is False

      --kvals                     Values of k for SPAdes assemblies. SPAdes needs to be 
                                  compiled to handle larger k-values! Default is 
                                  auto-detection by SPAdes

      --thresh                    Percent identity threshold for retaining Exonerate
                                  hits. Default is 55, but increase this if you are 
                                  worried about contaminant sequences

      --paralog_min_length_percentage <decimal> 
                                  Minimum length percentage of a SPAdes contig vs 
                                  reference protein query for a paralog warning to be 
                                  generated and a putative paralog contig to be 
                                  recovered. Default is 0.75 

      --depth_multiplier          Assign a long paralog as the "main" sequence if it 
                                  has a coverage depth <depth_multiplier> times all 
                                  other long paralogs. Set to zero to not use depth. 
                                  Default is 10

      --timeout_assemble          Kill long-running gene assemblies if they take longer 
                                  than X percent of average

      --timeout_exonerate_contigs Kill long-running processes if they take longer than 
                                  X seconds. Default is 120

      --target                    Use the target file sequence with this taxon name in 
                                  Exonerate searches for each gene. Other targets for 
                                  that gene will be used only for read sorting. Can be a 
                                  tab-delimited file (one <gene>\\t<taxon_name> per line) 
                                  or a single taxon name

      --exclude                   Do not use any sequence with the specified taxon name 
                                  string in Exonerate searches. Sequenced from this 
                                  taxon will still be used for read sorting

      --no_stitched_contig        Do not create stitched contigs; use longest Exonerate 
                                  hit only. Default is off



      --chimera_test_memory <int> Memory (RAM) amount in MB to use for bbmap.sh when
                                  peforming stitched-contig chimera tests. Default is 
                                  1000 MB

      --bbmap_subfilter <int>     Ban alignments with more than this many 
                                  substitutions when performing read-pair mapping to 
                                  supercontig reference (bbmap.sh). Default is 7

      --chimeric_stitched_contig_edit_distance <int>    
                                  Minimum number of base differences between one read 
                                  of a read pair vs the stitched-contig reference for a 
                                  read pair to be flagged as discordant. Default is 5
      
      --chimeric_stitched_contig_discordant_reads_cutoff <int>           
                                  Minimum number of discordant reads pairs required 
                                  to flag a stitched-contig as a potential chimera of 
                                  contigs from multiple paralogs. Default is 5

      --merged                    Merge forward and reverse reads, and run SPAdes 
                                  assembly with merged and unmerged (the latter 
                                  in interleaved format) data. Default is off
      
      --run_intronerate           Run the intronerate() fuction to recover intron 
                                  and supercontig sequences. Default is off, and so 
                                  fasta files in `subfolders 09_sequences_intron` and 
                                  `10_sequences_supercontig` will be empty

      --keep_intermediate_files   Keep all intermediate files and logs, which can be 
                                  useful for debugging. Default action is to delete 
                                  them, which greatly reduces the total file number

      --no_padding_supercontigs   If Intronerate is run, and a supercontig is created 
                                  by concatenating multiple SPAdes contigs, do not add 
                                  10 "N" characters between contig joins. By default, 
                                  Ns will be added

      --verbose_logging           If supplied, enable verbose login. NOTE: this can 
                                  increase the size of the log files by an order of 
                                  magnitude

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

/* 
Include a few default params here to print useful help (if requested) or if minimal input is not provided.
*/
params.help = false
params.illumina_reads_directory = false
params.target_file = false

// Check that input directories are provided
if (params.help || !params.illumina_reads_directory || (!params.targetfile_dna && !params.targetfile_aa)) {
  helpMessage()
  exit 0
}


// Check that paralog_warning_min_len_percent value is a decimal between 0 and 1
if (params.paralog_warning_min_len_percent < 0 || params.paralog_warning_min_len_percent >1) {
println("""
  The value for --paralog_warning_min_len_percent should be between 0 and 1. 
  Your value is ${params.paralog_warning_min_len_percent}""".stripIndent())
exit 0
}

// Check that non-overlapping options are provided
if (params.single_end && params.paired_and_single) {
  println('Please use --single_end OR --paired_and_single, not both!')
  exit 0
}
if (params.target_file_dna && params.target_file_aa) {
  println('Please use --target_file_dna OR --target_file_aa, not both!')
  exit 0
}
if (params.target_file_aa && params.use_bwa) {
  println('You can not use BWA with a target file containing protein sequences. \
  Please use BLASTx or DIAMOND, or provide a target file with nucleotide sequences.')
  exit 0
}


// Don't allow params.paired_and_single and params.use_trimmomatic
if (params.paired_and_single && params.use_trimmomatic) {
  println("""
    Trimmomatic can't be used with paired plus single reads yet - 
    let me know if this would be useful!""".stripIndent())
  exit 0
}



// Check for unrecognised pararmeters
allowed_params = ["no_stitched_contigs", "chimera_test_memory","chimeric_stitched_contig_edit_distance", \
"chimeric_stitched_contig_discordant_reads_cutoff", "merged", "paired_and_single", "single_end", "outdir", \
"illumina_reads_directory", "target_file", "help", "memory", "read_pairs_pattern", \
"single_pattern", "use_blastx", "num_forks", "cov_cutoff", "blastx_evalue", \
"paralog_warning_min_len_percent", "translate_target_file_for_blastx", "use_trimmomatic", \
"trimmomatic_leading_quality", "trimmomatic_trailing_quality", "trimmomatic_min_length", \
"trimmomatic_sliding_window_size", "trimmomatic_sliding_window_quality", "run_intronerate", \
"bbmap_subfilter", "combine_read_files", "combine_read_files_num_fields", "namelist", \
"keep_intermediate_files", "distribute_hi_mem", "use_diamond", "diamond_sensitivity", "single_cell_assembly",\
"timeout_assemble", "timeout_exonerate_contigs", "target", "exclude", "no_padding_supercontigs",\
"verbose_logging", "targetfile_aa", "targetfile_dna", "bwa"]

params.each { entry ->
  if (! allowed_params.contains(entry.key)) {
      println("The parameter <${entry.key}> is not known");
      exit 0;
  }
}


//////////////////////////////////
//  Target gene sequences file  //
//////////////////////////////////


if (params.targetfile_dna) {
  Channel
    .fromPath("${params.targetfile_dna}", checkIfExists: true)
    .first()
    .set { target_file_ch }
} else if (params.targetfile_aa) {
  Channel
    .fromPath("${params.targetfile_aa}", checkIfExists: true)
    .first()
    .set { target_file_ch }
}

end_field = params.combine_read_files_num_fields - 1  // Due to zero-based indexing

def getLibraryId( prefix ){
  /* 
  Function for grouping reads from multiple lanes, based on a shared filename 
  prefix preceeding the first underscore.
  */

  filename_list = prefix.split("_")
  groupby_select = filename_list[0..end_field]
  groupby_joined = groupby_select.join("_")
}


/////////////////////////////////////////////////////////
//  Create 'namelist.txt' file and associated channel  //
/////////////////////////////////////////////////////////

def user_provided_namelist_for_filtering = []

if (params.namelist) {
  user_provided_namelist_file = file("${params.namelist}", checkIfExists: true)
    .readLines()
    .each { user_provided_namelist_for_filtering << it }
  Channel
  .fromPath("${params.namelist}", checkIfExists: true)
  .first()
  .set { namelist_ch }

} else if (!params.single_end && !params.combine_read_files) {
  Channel
  .fromFilePairs("${params.illumina_reads_directory}/*_{$params.read_pairs_pattern}*.{fastq.gz,fastq,fq.gz,fq}", \
    flat : true, checkIfExists: true)
  .collectFile(name: "${params.outdir}/01_namelist/namelist.txt") { item -> item[0] + "\n" }
  .first()
  .set { namelist_ch }

} else if (!params.single_end && params.combine_read_files) {
  Channel
  .fromFilePairs("${params.illumina_reads_directory}/*_{$params.read_pairs_pattern}*.{fastq.gz,fastq,fq.gz,fq}", \
    flat : true, checkIfExists: true)
  .map { prefix, file1, file2 -> tuple(getLibraryId(prefix), file1, file2) }
  .groupTuple(sort:true)
  .collectFile(name: "${params.outdir}/01_namelist/namelist.txt") { item -> item[0] + "\n" }
  .first()
  .set { namelist_ch }

} else if (params.single_end && !params.combine_read_files) {
  Channel
  .fromPath("${params.illumina_reads_directory}/*_{$params.single_pattern}*.{fastq.gz,fastq,fq.gz,fq}", \
    checkIfExists: true)
  .map { file -> file.baseName.split("_${params.single_pattern}")[0] } // THIS NEEDS TO BE UNIQUE
  .unique()
  .collectFile(name: "${params.outdir}/01_namelist/namelist.txt", newLine: true)
  .first()
  .set { namelist_ch }

} else if (params.single_end && params.combine_read_files) {
  Channel
  .fromPath("${params.illumina_reads_directory}/*_{$params.single_pattern}*.{fastq.gz,fastq,fq.gz,fq}", \
    checkIfExists: true)
  .map { file -> tuple((file.baseName.split('_')[0..end_field]).join("_"), file) }
  .groupTuple(sort:true)
  .collectFile(name: "${params.outdir}/01_namelist/namelist.txt") { item -> item[0] + "\n" }
  .first()
  .set { namelist_ch }
}



if (user_provided_namelist_for_filtering) {
  user_provided_namelist_for_filtering = user_provided_namelist_for_filtering.findAll { item -> !item.isEmpty() }

  log.info("""
    INFO: A namelist has been supplied by the user. Only the following samples will be processed: ${user_provided_namelist_for_filtering}\n""".stripIndent())
}


//////////////////////////////
//  Illumina reads channel  //
//////////////////////////////

/*
Single-end reads.
Don't group reads from multi-lane (default).
*/
if (params.single_end && !params.combine_read_files && user_provided_namelist_for_filtering) {
  Channel
  .fromPath("${params.illumina_reads_directory}/*_{$params.single_pattern}*.{fastq.gz,fastq,fq.gz,fq}", \
    checkIfExists: true)
  .map { file -> tuple(file.baseName.split("_${params.single_pattern}")[0], file) } // THIS NEEDS TO BE UNIQUE
  .filter { it[0] in user_provided_namelist_for_filtering }
  // .view()
  .set { illumina_reads_single_end_ch }

} else if (params.single_end && !params.combine_read_files && 
  !user_provided_namelist_for_filtering) {
  Channel
  .fromPath("${params.illumina_reads_directory}/*_{$params.single_pattern}*.{fastq.gz,fastq,fq.gz,fq}", \
    checkIfExists: true)
  .map { file -> tuple(file.baseName.split("_${params.single_pattern}")[0], file) } // THIS NEEDS TO BE UNIQUE
  .set { illumina_reads_single_end_ch }

} else if (params.single_end && params.combine_read_files && user_provided_namelist_for_filtering) {
  Channel
  .fromPath("${params.illumina_reads_directory}/*_{$params.single_pattern}*.{fastq.gz,fastq,fq.gz,fq}", \
  checkIfExists: true)
  .map { file -> tuple((file.baseName.split('_')[0..end_field]).join("_"), file) }
  .groupTuple(sort:true)
  // .view()
  .filter { it[0] in user_provided_namelist_for_filtering }
  // .view()
  .set { illumina_reads_single_end_ch }

} else if (params.single_end && params.combine_read_files && !user_provided_namelist_for_filtering) {
  Channel
  .fromPath("${params.illumina_reads_directory}/*_{$params.single_pattern}*.{fastq.gz,fastq,fq.gz,fq}", \
  checkIfExists: true)
  .map { file -> tuple((file.baseName.split('_')[0..end_field]).join("_"), file) }
  .groupTuple(sort:true)
  .set { illumina_reads_single_end_ch }

} else {
  illumina_reads_single_end_ch = Channel.empty()
}


/*
Paired-end reads and a file of unpaired reads.
*/
if (params.paired_and_single) {
  Channel
  .fromFilePairs("${params.illumina_reads_directory}/*_{$params.read_pairs_pattern,$params.single_pattern}*.{fastq.gz,fastq,fq.gz,fq}", flat : true,
  checkIfExists: true, size: 3)
  .set { illumina_paired_reads_with_unpaired_ch }
} else {
  illumina_paired_reads_with_unpaired_ch = Channel.empty()
}


/* 
Paired-end reads.
Don't group reads from multi-lane (default).
*/
if (!params.paired_and_single && !params.single_end  && !params.combine_read_files && user_provided_namelist_for_filtering) {
  Channel
    .fromFilePairs("${params.illumina_reads_directory}/*_{$params.read_pairs_pattern}*.{fastq.gz,fastq,fq.gz,fq}", \
    flat : true, checkIfExists: true)
    .filter { it[0] in user_provided_namelist_for_filtering }
    // .view()
    .set { illumina_paired_reads_ch }

} else if (!params.paired_and_single && !params.single_end  && !params.combine_read_files && !user_provided_namelist_for_filtering) {
    Channel
    .fromFilePairs("${params.illumina_reads_directory}/*_{$params.read_pairs_pattern}*.{fastq.gz,fastq,fq.gz,fq}", \
    flat : true, checkIfExists: true)
    // .view()
    .set { illumina_paired_reads_ch }

} else if (!params.paired_and_single && !params.single_end  && params.combine_read_files && user_provided_namelist_for_filtering) {
    Channel
    .fromFilePairs("${params.illumina_reads_directory}/*_{$params.read_pairs_pattern}*.{fastq.gz,fastq,fq.gz,fq}", \
    flat : true, checkIfExists: true)
    .map { prefix, file1, file2 -> tuple(getLibraryId(prefix), file1, file2) }
    .groupTuple(sort:true)
    .filter { it[0] in user_provided_namelist_for_filtering }
    // .view()
    .set { illumina_paired_reads_ch }

} else if (!params.paired_and_single && !params.single_end && params.combine_read_files && !user_provided_namelist_for_filtering) {
    Channel
    .fromFilePairs("${params.illumina_reads_directory}/*_{$params.read_pairs_pattern}*.{fastq.gz,fastq,fq.gz,fq}", \
    flat : true, checkIfExists: true)
    .map { prefix, file1, file2 -> tuple(getLibraryId(prefix), file1, file2) }
    .groupTuple(sort:true)
    .set { illumina_paired_reads_ch }

} else {
  illumina_paired_reads_ch = Channel.empty()
}


/////////////////////////////
//  DEFINE DSL2 PROCESSES  //
/////////////////////////////

process COMBINE_LANES_PAIRED_END {
  /*
  If `--combine_read_files` flag is set, combine lanes when using paired-end R1 and R2 reads.
  */

  label 'in_container'
  // echo true
  publishDir "$params.outdir/02_reads_combined_lanes", mode: 'copy', pattern: "*.fastq*"

  if (params.num_forks) {
    maxForks params.num_forks
  }

  when:
    params.combine_read_files

  input:
    tuple val(prefix), path(reads_R1), path(reads_R2)

  output:
    tuple val(prefix), path("*R1.fastq*"), path("*R2.fastq*"), emit: combined_lane_paired_reads

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


process COMBINE_LANES_SINGLE_END {
  /*
  If `--combine_read_files` flag is set, combine lanes when using single-end reads only,
  */

  label 'in_container'
  // echo true
  publishDir "$params.outdir/02_reads_combined_lanes", mode: 'copy', pattern: "*.fastq*"

  if (params.num_forks) {
    maxForks params.num_forks
  }

  when:
    params.combine_read_files

  input:
    tuple val(prefix), path(reads_single)

  output:
    tuple val(prefix), path("*single.fastq*"), emit: combined_lane_single_reads_ch

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


process TRIMMOMATIC_PAIRED {
  /*
  If `--use_trimmomatic` flag is set, run optional Trimmomatic step for paired-end reads.
  */

  // echo true
  label 'in_container'
  publishDir "$params.outdir/03a_trimmomatic_logs", mode: 'copy', pattern: "*.log"
  publishDir "$params.outdir/03b_trimmomatic_paired_and_single_reads", mode: 'copy', pattern: "*_paired.fq*"
  publishDir "$params.outdir/03b_trimmomatic_paired_and_single_reads", mode: 'copy', pattern: "*_R1-R2_unpaired.fq*"


  if (params.num_forks) {
    maxForks params.num_forks
  }

  when:
    params.use_trimmomatic

  input:
    tuple val(prefix), path(reads_R1), path(reads_R2)

  output:
    path("*")
    tuple val(prefix), path("*R1_paired*"), path("*R2_paired*"), path("*R1-R2_unpaired*"), emit: trimmed_paired_and_orphaned_ch

  script:
    read_pairs_pattern_list = params.read_pairs_pattern?.tokenize(',')

    """
    R1=${reads_R1}
    R2=${reads_R2}
    sampleID_R1=\${R1%_${read_pairs_pattern_list[0]}*}
    sampleID_R2=\${R2%_${read_pairs_pattern_list[1]}*}

    echo \$R1
    echo \$R2

    if [[ \$R1 = *.gz ]]
      then 
        R1_filename_strip_gz="\${R1%.gz}"
        fastq_extension="\${R1_filename_strip_gz##*.}"

        output_forward_paired=\${sampleID_R1}_R1_paired.fq.gz
        output_reverse_paired=\${sampleID_R2}_R2_paired.fq.gz
        output_forward_unpaired=\${sampleID_R1}_R1_unpaired.fq.gz
        output_reverse_unpaired=\${sampleID_R2}_R2_unpaired.fq.gz
        output_both_unpaired=\${sampleID_R1}_R1-R2_unpaired.fq.gz

      else
        fastq_extension="\${R1##*.}"

        output_forward_paired=\${sampleID_R1}_R1_paired.fq
        output_reverse_paired=\${sampleID_R2}_R2_paired.fq
        output_forward_unpaired=\${sampleID_R1}_R1_unpaired.fq
        output_reverse_unpaired=\${sampleID_R2}_R2_unpaired.fq
        output_both_unpaired=\${sampleID_R1}_R1-R2_unpaired.fq
    fi

    # Write adapters fasta file:
    echo -e ">PrefixPE/1\nTACACTCTTTCCCTACACGACGCTCTTCCGATCT\n>PrefixPE/2\nGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT\n>PE1\nTACACTCTTTCCCTACACGACGCTCTTCCGATCT\n>PE1_rc\nAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA\n>PE2\nGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT\n>PE2_rc\nAGATCGGAAGAGCACACGTCTGAACTCCAGTCA" > TruSeq3-PE-2.fa

    # Run Trimmomtic:
    trimmomatic PE -phred33 -threads ${task.cpus} \
    ${reads_R1} ${reads_R2} \${output_forward_paired} \${output_forward_unpaired} \
    \${output_reverse_paired} \${output_reverse_unpaired} \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:1:true \
    LEADING:${params.trimmomatic_leading_quality} \
    TRAILING:${params.trimmomatic_trailing_quality} \
    SLIDINGWINDOW:${params.trimmomatic_sliding_window_size}:${params.trimmomatic_sliding_window_quality} \
    MINLEN:${params.trimmomatic_min_length} 2>&1 | tee \${sampleID_R1}.log 
    cat \${output_forward_unpaired} \${output_reverse_unpaired} > \${output_both_unpaired}
    """
}


process TRIMMOMATIC_SINGLE {
  /*
  If `--use_trimmomatic` flag is set, run optional Trimmomatic step for single-end reads.
  */

  // echo true
  label 'in_container'
  publishDir "$params.outdir/03a_trimmomatic_logs", mode: 'copy', pattern: "*.log"
  publishDir "$params.outdir/03c_trimmomatic_single_reads", mode: 'copy', pattern: "*_single*"

  if (params.num_forks) {
    maxForks params.num_forks
  }

  when:
    params.use_trimmomatic

  input:
    tuple val(prefix), path(reads_single)

  output:
    path("*")
    tuple val(prefix), file("*single*"), emit: trimmed_single_ch

  script:
    """
    single=${reads_single}


    if [[ \$single = *.gz ]]
      then 
        output_single=${prefix}_trimmed_single.fq.gz
      else
        output_single=${prefix}_trimmed_single.fq
    fi

    echo -e ">TruSeq3_IndexedAdapter\nAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC\n>TruSeq3_UniversalAdapter\nAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA\n" > TruSeq3-SE.fa
    trimmomatic SE -phred33 -threads ${task.cpus} \
    ${reads_single} \${output_single} ILLUMINACLIP:TruSeq3-SE.fa:2:30:10:1:true \
    LEADING:${params.trimmomatic_leading_quality} \
    TRAILING:${params.trimmomatic_trailing_quality} \
    SLIDINGWINDOW:${params.trimmomatic_sliding_window_size}:${params.trimmomatic_sliding_window_quality} \
    MINLEN:${params.trimmomatic_min_length} 2>&1 | tee ${prefix}.log 
    """
}



process ASSEMBLE_SINGLE_END {
  /*
  Run the `hybpiper assemble` command for input files: [single_end]
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}/06_summary_stats", mode: 'copy', pattern: "${prefix}/${prefix}_genes_with_supercontigs.csv"
  publishDir "${params.outdir}/06_summary_stats", mode: 'copy', pattern: "${prefix}/${prefix}_supercontigs_with_discordant_reads.csv"

  if (params.num_forks) {
    maxForks params.num_forks
  }

  when:
    params.single_end

  input:
    path(target_file)
    tuple val(prefix), path(reads_single)

  output:
    path("${prefix}"), emit: assemble_with_single_end_ch optional true
    path("${prefix}/${prefix}_genes_with_supercontigs.csv") optional true
    path("${prefix}/${prefix}_supercontigs_with_discordant_reads.csv") optional true

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
      cleanup = "python /HybPiper/cleanup.py ${prefix}"
    } else {
      cleanup = ''
    }
    assemble_command = "python /HybPiper/assemble.py -b ${target_file} -r ${reads_single} --prefix ${prefix} --cpu ${task.cpus} " + command_list.join(' ')

    """
    echo ${assemble_command}
    ${assemble_command}
    ${cleanup}
    """
}


process ASSEMBLE_PAIRED_AND_SINGLE_END {
  /*
  Run the `hybpiper assemble` command for input files: [R1, R1, R1-R2_unpaired]
  */

  //echo true
  label 'in_container'
  publishDir "${params.outdir}/06_summary_stats", mode: 'copy', pattern: "${pair_id}/${pair_id}_genes_with_stitched_contig.csv"
  publishDir "${params.outdir}/06_summary_stats", mode: 'copy', pattern: "${pair_id}/${pair_id}_genes_derived_from_putative_chimeric_stitched_contig.csv"
  publishDir "${params.outdir}/06_summary_stats", mode: 'copy', pattern: "${pair_id}/${pair_id}_genes_with_long_paralog_warnings.txt"
  publishDir "${params.outdir}/06_summary_stats", mode: 'copy', pattern: "${pair_id}/${pair_id}_genes_with_paralog_warnings_by_contig_depth.csv"

  if (params.num_forks) {
    maxForks params.num_forks
  }

  when:
    (params.use_trimmomatic || params.paired_and_single)

  input:
    path(target_file) 
    tuple val(pair_id), path(reads_R1), path(reads_R2), path(reads_unpaired)

  output:
    path("${pair_id}"), emit: assemble_with_unPaired_ch optional true
    path("${pair_id}/${pair_id}_genes_with_stitched_contig.csv") optional true
    path("${pair_id}/${pair_id}_genes_derived_from_putative_chimeric_stitched_contig.csv") optional true
    path("${pair_id}/${pair_id}_genes_with_long_paralog_warnings.txt") optional true
    path("${pair_id}/${pair_id}_genes_with_paralog_warnings_by_contig_depth.csv") optional true

  script:
    def command_list = []

    if (params.nosupercontigs) {
      command_list << "--nosupercontigs"
      }
    if (params.memory) {
      command_list << "--memory ${params.memory}"
      }
    if (params.bbmap_subfilter) {
      command_list << "--bbmap_subfilter ${params.bbmap_subfilter}"
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
      cleanup = "python /HybPiper/cleanup.py ${pair_id}"
    } else {
      cleanup = ''
    }
    assemble_command = "python /HybPiper/assemble.py -b ${target_file} -r ${reads_R1} ${reads_R2} --unpaired ${reads_unpaired} --prefix ${pair_id} --cpu ${task.cpus} " + command_list.join(' ')

    script:
    """
    echo ${assemble_command}
    ${assemble_command}
    ${cleanup}
    """
  } 


process ASSEMBLE_PAIRED_END {
  /*
  Run the `hyvpiper assemble` command for input files: [R1, R1]
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}/06_summary_stats", mode: 'copy', pattern: "${pair_id}/${pair_id}_genes_with_stitched_contig.csv"
  publishDir "${params.outdir}/06_summary_stats", mode: 'copy', pattern: "${pair_id}/${pair_id}_genes_derived_from_putative_chimeric_stitched_contig.csv"
  publishDir "${params.outdir}/06_summary_stats", mode: 'copy', pattern: "${pair_id}/${pair_id}_genes_with_long_paralog_warnings.txt"
  publishDir "${params.outdir}/06_summary_stats", mode: 'copy', pattern: "${pair_id}/${pair_id}_genes_with_paralog_warnings_by_contig_depth.csv"

  if (params.num_forks) {
    maxForks params.num_forks
  }

  when:
    (!params.paired_and_single && !params.single_end && !params.use_trimmomatic)

  input:
    path(target_file) 
    tuple val(pair_id), path(reads_R1), path(reads_R2)

  output:
    path("${pair_id}"), emit: assemble_ch optional true
    path("${pair_id}/${pair_id}_genes_with_stitched_contig.csv") optional true
    path("${pair_id}/${pair_id}_genes_derived_from_putative_chimeric_stitched_contig.csv") optional true
    path("${pair_id}/${pair_id}_genes_with_long_paralog_warnings.txt") optional true
    path("${pair_id}/${pair_id}_genes_with_paralog_warnings_by_contig_depth.csv") optional true

  script:
    def command_list = []

    if (params.targetfile_dna) {
      command_list << "--targetfile_dna ${params.targetfile_dna}"
      }
    if (params.targetfile_aa) {
      command_list << "--targetfile_dna ${params.targetfile_dna}"
      }
    if (params.bwa) {
      command_list << "--bwa"
      }
    if (params.use_diamond) {
      command_list << "--diamond"
      }
    if (params.diamond_sensitivity) {
      command_list << "--diamond_sensitivity ${params.diamond_sensitivity}"
      }
    if (params.distribute_hi_mem) {
      command_list << "--distribute_hi_mem"
      }
    if (params.evalue) {
      command_list << "--evalue ${params.evalue}"
      }
    if (params.max_target_seqs) {
      command_list << "--max_target_seqs ${params.max_target_seqs}"
      }
    if (params.cov_cutoff) {
      command_list << "--cov_cutoff ${params.cov_cutoff}"
      }
    if (params.single_cell_assembly) {
      command_list << "--single_cell_assembly"
      }
    if (params.kvals) {
      command_list << "--kvals ${params.kvals}"
      }
    if (params.paralog_min_length_percentage) {
      command_list << "--paralog_min_length_percentage ${params.paralog_min_length_percentage}"
      }
    if (params.timeout_assemble) {
      command_list << "--timeout_assemble ${params.timeout_assemble}"
      }
    if (params.timeout_exonerate_contigs) {
      command_list << "--timeout_exonerate_contigs ${params.timeout_exonerate_contigs}"
      }
    if (params.target) {
      command_list << "--target ${params.target}"
      }
    if (params.exclude) {
      command_list << "--exclude ${params.exclude}"
      }
    if (params.no_stitched_contig) {
      command_list << "--no_stitched_contig"
      }
    if (params.chimera_test_memory) {
      command_list << "--bbmap_memory ${params.chimera_test_memory}"
      }
    if (params.bbmap_subfilter) {
      command_list << "--bbmap_subfilter ${params.bbmap_subfilter}"
      }
    if (params.chimeric_stitched_contig_edit_distance) {
      command_list << "--chimeric_stitched_contig_edit_distance ${params.chimeric_stitched_contig_edit_distance}"
      }
    if (params.chimeric_stitched_contig_discordant_reads_cutoff) {
      command_list << "--chimeric_stitched_contig_discordant_reads_cutoff ${params.chimeric_stitched_contig_discordant_reads_cutoff}"
      }
    if (params.merged) {
      command_list << "--merged"
      }
    if (params.run_intronerate) {
      command_list << "--run_intronerate"
      }
    if (params.keep_intermediate_files) {
      command_list << "--keep_intermediate_files"
      }
    if (params.no_padding_supercontigs) {
      command_list << "--no_padding_supercontigs"
      }
    if (params.verbose_logging) {
      command_list << "--verbose_logging"
      }

    assemble_command = "hybpiper assemble -r ${reads_R1} ${reads_R2} --prefix ${pair_id} --cpu ${task.cpus} " + command_list.join(' ')


    script:
    """
    echo "Executing command: ${assemble_command}"
    ${assemble_command}
    """
}


process SUMMARY_STATS {
/*
Run the `hybpiper stats` command.
*/

  // echo true
  label 'in_container'
  publishDir "${params.outdir}/06_summary_stats", mode: 'copy'

  input:
    path(assemble)
    path(target_file)
    path(namelist)

  output:
    path("hybpiper_stats.tsv"), emit: stats_file
    path("seq_lengths.tsv"), emit: seq_lengths_file

  script:
  if (params.targetfile_dna) {
  """
  hybpiper stats -t_dna ${target_file} gene ${namelist}
  """
  } else if (params.targetfile_dna) {
  """
  hybpiper stats -t_aa ${target_file} gene ${namelist}
  """
  }

}


process VISUALISE {
  /*
  Run the `hybpiper recovery_heatmap` command.
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}/05_visualise", mode: 'copy'

  input:
    path(seq_lengths_file)

  output:
    path("recovery_heatmap.png")

  script:
    """
    hybpiper recovery_heatmap ${seq_lengths_file}
    """
}


process RETRIEVE_SEQUENCES {
  /*
  Run the hybpiper retrieve_sequences command for all sequence types.
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}/07_sequences_dna", mode: 'copy', pattern: "*.FNA"
  publishDir "${params.outdir}/08_sequences_aa", mode: 'copy', pattern: "*.FAA"
  publishDir "${params.outdir}/09_sequences_intron", mode: 'copy', pattern: "*introns.fasta"
  publishDir "${params.outdir}/10_sequences_supercontig", mode: 'copy', pattern: "*supercontig.fasta"

  input:
    path(assemble)
    path(target_file)
    path(namelist)

  output:
    path("*.FNA")
    path("*.FAA")
    path("*.fasta")

  script:
  if (params.targetfile_dna) {
    """
    hybpiper retrieve_sequences -t_dna ${target_file} --sample_names ${namelist} dna
    hybpiper retrieve_sequences -t_dna ${target_file} --sample_names ${namelist} aa
    hybpiper retrieve_sequences -t_dna ${target_file} --sample_names ${namelist} intron
    hybpiper retrieve_sequences -t_dna ${target_file} --sample_names ${namelist} supercontig
    """
  } else if (params.targetfile_aa) {
    """
    hybpiper retrieve_sequences -t_aa ${target_file} --sample_names ${namelist} . dna
    hybpiper retrieve_sequences -t_aa ${target_file} --sample_names ${namelist} . aa
    hybpiper retrieve_sequences -t_aa ${target_file} --sample_names ${namelist} . intron
    hybpiper retrieve_sequences -t_aa ${target_file} --sample_names ${namelist} . supercontig
    """
  }
}


process PARALOG_RETRIEVER {
  /*
  Run hybpiper `paralog_retriever` command.
  */

  //echo true
  label 'in_container'
  publishDir "${params.outdir}/05_visualise", mode: 'copy', pattern: "paralog_heatmap.png"
  publishDir "${params.outdir}/11_paralogs", mode: 'copy', pattern: "paralogs_all/*paralogs_all.fasta", saveAs: { filename -> file(filename).getName() }
  publishDir "${params.outdir}/12_paralogs_no_chimeras", mode: 'copy', pattern: "paralogs_no_chimeras/*paralogs_no_chimeras.fasta", saveAs: { filename -> file(filename).getName() }
  publishDir "${params.outdir}/11_paralogs/logs", mode: 'copy', pattern: "*report*"

  input:
    // path(paralog_complete_list)
    path(assemble)
    path(namelist)
    path(target_file)


  output:
    path("paralogs_all/*paralogs_all.fasta")
    path("paralogs_no_chimeras/*paralogs_no_chimeras.fasta")
    path("paralog_report.tsv")
    path("paralogs_above_threshold_report.txt") 

  script:
    if (params.targetfile_dna) {
      """
      hybpiper paralog_retriever  ${namelist} -t_dna ${target_file}
      """
    } else if (params.targetfile_aa) {
      """
      hybpiper paralog_retriever  ${namelist} -t_aa ${target_file}
      """
    }
}



////////////////////////
//  Define workflows  //
////////////////////////

workflow {

  // // Run OPTIONAL translate target file step:
  // TRANSLATE_TARGET_FILE( target_file_ch )

  // // Set up input channel for target file:
  // if (!params.translate_target_file_for_blastx) {
  //   target_file_ch = target_file_ch
  // } else {
  //   target_file_ch = TRANSLATE_TARGET_FILE.out.translated_target_file
  // }

  // Run OPTIONAL combine read file step: 
  COMBINE_LANES_PAIRED_END( illumina_paired_reads_ch )
  COMBINE_LANES_SINGLE_END( illumina_reads_single_end_ch )

  // Set up correct channel for combined vs non-combined:
  if (params.combine_read_files) {
    trimmomatic_PE_input_ch = COMBINE_LANES_PAIRED_END.out.combined_lane_paired_reads
    trimmomatic_SE_input_ch = COMBINE_LANES_SINGLE_END.out.combined_lane_single_reads_ch
  } else {
    trimmomatic_PE_input_ch = illumina_paired_reads_ch
    trimmomatic_SE_input_ch = illumina_reads_single_end_ch
  }

  // Run OPTIONAL trimmomatic QC step:
  TRIMMOMATIC_PAIRED( trimmomatic_PE_input_ch )
  TRIMMOMATIC_SINGLE( trimmomatic_SE_input_ch )

  // Set up input channels for assemble.py:
  if (params.use_trimmomatic) {
    assemble_with_single_end_only_input_ch = TRIMMOMATIC_SINGLE.out.trimmed_single_ch
    assemble_with_unpaired_input_ch = TRIMMOMATIC_PAIRED.out.trimmed_paired_and_orphaned_ch
    assemble_no_unpaired_input_ch = Channel.empty()
  } else if (params.combine_read_files) {
    assemble_with_single_end_only_input_ch = COMBINE_LANES_SINGLE_END.out.combined_lane_single_reads_ch
    assemble_with_unpaired_input_ch = Channel.empty()
    assemble_no_unpaired_input_ch = COMBINE_LANES_PAIRED_END.out.combined_lane_paired_reads
  } else {
    assemble_with_single_end_only_input_ch = illumina_reads_single_end_ch
    assemble_with_unpaired_input_ch = illumina_paired_reads_with_unpaired_ch
    assemble_no_unpaired_input_ch = illumina_paired_reads_ch
  }

  // Run hybpiper assemble:
  ASSEMBLE_PAIRED_AND_SINGLE_END( target_file_ch, assemble_with_unpaired_input_ch )
  ASSEMBLE_PAIRED_END( target_file_ch, assemble_no_unpaired_input_ch )
  ASSEMBLE_SINGLE_END ( target_file_ch, assemble_with_single_end_only_input_ch )

  // Run hybpiper stats:
  SUMMARY_STATS( ASSEMBLE_PAIRED_AND_SINGLE_END.out.assemble_with_unPaired_ch.collect().mix(ASSEMBLE_PAIRED_END.out.assemble_ch).collect().mix(ASSEMBLE_SINGLE_END.out.assemble_with_single_end_ch).collect(), target_file_ch, namelist_ch ) 

  // Run hybpiper recovery_heatmap: 
  VISUALISE( SUMMARY_STATS.out.seq_lengths_file ) 

  // Run retrieve_sequences.py script for all sequence types:
  RETRIEVE_SEQUENCES( ASSEMBLE_PAIRED_AND_SINGLE_END.out.assemble_with_unPaired_ch.collect().mix(ASSEMBLE_PAIRED_END.out.assemble_ch).collect().mix(ASSEMBLE_SINGLE_END.out.assemble_with_single_end_ch).collect(), target_file_ch, namelist_ch )

  // Run hybpiper paralog_retriever: 
  PARALOG_RETRIEVER( ASSEMBLE_PAIRED_AND_SINGLE_END.out.assemble_with_unPaired_ch.collect().mix(ASSEMBLE_PAIRED_END.out.assemble_ch).collect().mix(ASSEMBLE_SINGLE_END.out.assemble_with_single_end_ch).collect(), namelist_ch, target_file_ch )
} 

///////////////////////////////////////////////////
/////////////////  End of script  /////////////////
///////////////////////////////////////////////////
