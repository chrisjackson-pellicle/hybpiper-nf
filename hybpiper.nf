#!/usr/bin/env nextflow

/////////////////////////////////////////////////////////////////////////////////////////
////////////////////  Nextflow Pipeline for HybPiper version 2.3.1  /////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

nextflow.enable.dsl=2

/////////////////////////////////////////////////////////////////////////////////////////
// PRINT HYBPIPER-NF NEXTFLOW SCRIPT VERSION
/////////////////////////////////////////////////////////////////////////////////////////

if( params.remove('version') ) {
    println('hybpiper-nf version 1.1.1, running HybPiper version 2.3.2')
    exit 0
} 


/////////////////////////////////////////////////////////////////////////////////////////
// HELP MESSAGE IF DEFAULT WORKFLOW RUN (i.e. no entry point provided)
/////////////////////////////////////////////////////////////////////////////////////////

workflow {

  println(
  """
  Please provide a workflow name using the parameter "-entry" and one of:

      check_targetfile
      fix_targetfile
      assemble

  Use e.g. `-entry check_targetfile --help` to see the full parameters for the workflow.
  """
  )
  exit(0)
}


/////////////////////////////////////////////////////////////////////////////////////////
// MESSAGE IF WORKFLOW filter_by_length USED
/////////////////////////////////////////////////////////////////////////////////////////

if (workflow.commandLine.contains('-entry filter_by_length')) {
  println("The 'filter_by_length' subcommand is not supported yet.");
  exit 0
} 

/////////////////////////////////////////////////////////////////////////////////////////
// CHECK FOR UNRECOGNISED PARAMETERS
/////////////////////////////////////////////////////////////////////////////////////////

allowed_params = [
                  ////////////////////////////////////////////////////////////////////////
                  // Input:
                  ////////////////////////////////////////////////////////////////////////
                  "illumina_reads_directory",
                  "targetfile_aa", 
                  "targetfile_dna",
                  "read_pairs_pattern",
                  "single_pattern",
                  "paired_and_single", 
                  "single_end",
                  "combine_read_files", 
                  "combine_read_files_num_fields",

                  ////////////////////////////////////////////////////////////////////////
                  // Trimming reads (optional):
                  ////////////////////////////////////////////////////////////////////////
                  "use_trimmomatic",
                  "trimmomatic_leading_quality", 
                  "trimmomatic_trailing_quality",
                  "trimmomatic_min_length",
                  "trimmomatic_sliding_window_size", 
                  "trimmomatic_sliding_window_quality",

                  ////////////////////////////////////////////////////////////////////////
                  // `hybpiper assemble`
                  ////////////////////////////////////////////////////////////////////////

                  // Options for step: map_reads:
                  "bwa",
                  "diamond", 
                  "diamond_sensitivity",
                  "evalue",
                  "max_target_seqs",

                  // Options for step: distribute_reads:
                  "distribute_low_mem",

                  // Options for step: assemble_reads:
                  "cov_cutoff",
                  "single_cell_assembly",
                  "kvals",
                  "merged",
                  "timeout_assemble_reads",
                  // "no_spades_eta", NOT USED - specified as default below

                  // Options for step: extract_contigs:
                  "not_protein_coding",
                  "extract_contigs_blast_task",
                  "extract_contigs_blast_evalue",
                  "extract_contigs_blast_word_size",
                  "extract_contigs_blast_gapopen",
                  "extract_contigs_blast_gapextend",
                  "extract_contigs_blast_penalty",
                  "extract_contigs_blast_reward",
                  "extract_contigs_blast_perc_identity",
                  "extract_contigs_blast_max_target_seqs",
                  "thresh",
                  "paralog_min_length_percentage",
                  "depth_multiplier",
                  "target",
                  "exclude",
                  "timeout_extract_contigs",
                  "no_stitched_contig",
                  "no_pad_stitched_contig_gaps_with_n",
                  "chimeric_stitched_contig_check",
                  "bbmap_memory",
                  "bbmap_subfilter",
                  "bbmap_threads",
                  "chimeric_stitched_contig_edit_distance",
                  "chimeric_stitched_contig_discordant_reads_cutoff",
                  "trim_hit_sliding_window_size",
                  "trim_hit_sliding_window_thresh",
                  "exonerate_skip_hits_with_frameshifts",
                  "exonerate_skip_hits_with_internal_stop_codons",
                  "exonerate_skip_hits_with_terminal_stop_codons",
                  "exonerate_refine_full",
                  "no_intronerate",
                  "no_padding_supercontigs",

                  // General `hybpiper assemble` pipeline options:
                  // "prefix", NOT USED
                  // "start_from", NOT USED
                  // "end_with", NOT USED
                  // "force_overwrite", NOT USED
                  // "cpu", NOT USED - manually specified as task.cpus below
                  "compress_sample_folder",
                  "skip_targetfile_checks",
                  "keep_intermediate_files",
                  "verbose_logging",
                  // "hybpiper_output", NOT USED
                  "run_profiler",

                  // General hybpiper.nf pipeline options:
                  "namelist",
                  "outdir",
                  "help", 
                  "memory", 
                  "num_forks",
                  "chimera_test_memory",
                   
                  ////////////////////////////////////////////////////////////////////////
                  // `hybpiper check_targetfile`
                  ////////////////////////////////////////////////////////////////////////
                  "check_targetfile",
                  "no_terminal_stop_codons",
                  "sliding_window_size",
                  "complexity_minimum_threshold",

                  ////////////////////////////////////////////////////////////////////////
                  // `hybpiper fix_targetfile`
                  ////////////////////////////////////////////////////////////////////////
                  "control_file",
                  "allow_gene_removal", 
                  "reference_protein_file", 
                  "maximum_distance",  
                  "filter_by_length_percentage",  
                  "keep_low_complexity_sequences",
                  "alignments",  
                  "concurrent_alignments",  // CHECK IF USED!!  
                  "threads_per_concurrent_alignment", // CHECK IF USED!! 
                  "write_all_fasta_files",
                    
                  
                  "figure_length", 
                  "figure_height", 
                  "sample_text_size", 
                  "gene_text_size", 
                  "heatmap_filetype", 
                  "heatmap_dpi", 
                  "paralogs_list_threshold_percentage"
                  ]

params.each { entry ->
  if (!allowed_params.contains(entry.key)) {
      println("The parameter <${entry.key}> is not known");
      exit 0;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////
// SET CHANNELS FOR TARGET FILE AND CONTROL FILE
/////////////////////////////////////////////////////////////////////////////////////////

//  Set the target file channel:
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


//  Set the control *.ctl file channel:
if (params.control_file) {
  Channel
    .fromPath("${params.control_file}", checkIfExists: true)
    .first()
    .set { control_file_ch }
} 


/////////////////////////////////////////////////////////////////////////////////////////
// DEFINE SOME FUNCTIONS     
/////////////////////////////////////////////////////////////////////////////////////////

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


/////////////////////////////////////////////////////////////////////////////////////////
// CHECK TARGET FILE
/////////////////////////////////////////////////////////////////////////////////////////

def check_targetfile_help() {
  log.info """
  #######################################################################################
  #################### hybpiper check_targetfile options: ###############################
  #######################################################################################

      --targetfile_dna TARGETFILE_DNA, -t_dna TARGETFILE_DNA
                            FASTA file containing DNA target sequences for each
                            gene. The fasta headers must follow the naming
                            convention: >TaxonID-geneName
      --targetfile_aa TARGETFILE_AA, -t_aa TARGETFILE_AA
                            FASTA file containing amino-acid target sequences for
                            each gene. The fasta headers must follow the naming
                            convention: >TaxonID-geneName
      --no_terminal_stop_codons
                            When testing for open reading frames, do not allow a
                            translated frame to have a single stop codon at the
                            C-terminus of the translated protein sequence. Default
                            is False.
      --sliding_window_size SLIDING_WINDOW_SIZE
                            Number of characters (single-letter DNA or amino-acid
                            codes) to include in the sliding window when checking
                            for sequences with low-complexity-regions.
      --complexity_minimum_threshold COMPLEXITY_MINIMUM_THRESHOLD
                            Minimum threshold value. Beneath this value, the
                            sequence in the sliding window is flagged as low
                            complexity, and the corresponding target file sequence
                            is reported as having low-complexity regions.
      --run_profiler        If supplied, run the subcommand using cProfile. Saves
                            a *.csv file of results

  #######################################################################################
  """.stripIndent()
}

// Check for `--help` parameter if running `-entry check_targetfile`:
if (workflow.commandLine.contains('-entry check_targetfile') && params.remove('help') ) {
  check_targetfile_help()
  exit 0
} 

// Check minimal input provided if `-entry check_targetfile`:
if (workflow.commandLine.contains('-entry check_targetfile')  && 
(!params.targetfile_dna && !params.targetfile_aa)) {
  println(
  """
  ERROR: Parameter "-entry check_targetfile" provided, but no target file provided!
  Please provide your target file using the "--targetfile_dna" or "--targetfile_aa" parameter.
  """
  )
  exit 0
} 

// Create `hybpiper check_targetfile` command string:
def check_targetfile_command_list = []

if (params.targetfile_dna) {
  File target_file = new File(params.targetfile_dna)
  target_file_basename = target_file.getName()
  check_targetfile_command_list << "--targetfile_dna ${target_file_basename}"
  }
if (params.targetfile_aa) {
  File target_file = new File(params.targetfile_aa)
  target_file_basename = target_file.getName()
  check_targetfile_command_list << "--targetfile_aa ${target_file_basename}"
  }
if (params.sliding_window_size) {
  check_targetfile_command_list << "--sliding_window_size ${params.sliding_window_size}"
  } 
if (params.complexity_minimum_threshold) {
  check_targetfile_command_list << "--complexity_minimum_threshold ${params.complexity_minimum_threshold}"
  }
if (params.run_profiler) {
  check_targetfile_command_list << "--run_profiler ${params.run_profiler}"
  } 

// Workflow to run the check_targetfile_main workflow with target_file_ch as input:
workflow check_targetfile {
    check_targetfile_main( target_file_ch )
}

// Workflow to run the CHECK_TARGETFILE process
workflow check_targetfile_main {
    take: target_file
    main:
        CHECK_TARGETFILE(target_file)
        CHECK_TARGETFILE.out.check_results.view()
}


/////////////////////////////////////////////////////////////////////////////////////////
// FIX TARGET FILE
/////////////////////////////////////////////////////////////////////////////////////////

def fix_targetfile_help() {
  log.info """
  #######################################################################################
  ###################### hybpiper fix_targetfile options: ###############################
  #######################################################################################

      --control_file CONTROL_FILE
                            The *.ctl file, as output by the command 
                            "-entry check_targetfile".

      --targetfile_dna TARGETFILE_DNA, -t_dna TARGETFILE_DNA
                            FASTA file containing DNA target sequences for each gene. 
                            The fasta headers must follow the naming convention: 
                            >TaxonID-geneName
      --targetfile_aa TARGETFILE_AA, -t_aa TARGETFILE_AA
                            FASTA file containing amino-acid target sequences for each 
                            gene. The fasta headers must follow the naming convention: 
                            >TaxonID-geneName
      --no_terminal_stop_codons
                            When testing for open reading frames, do not allow a 
                            translated frame to have a single stop codon at the 
                            C-terminus of the translated protein sequence. Default is 
                            False. If supplied, this parameter will override the setting 
                            in the *.ctl file.
      --allow_gene_removal  Allow frame-correction and filtering steps to remove all 
                            representative sequences for a given gene. Default is False; 
                            HybPiper will exit with an information message instead. 
                            If supplied, this parameter will override the setting in 
                            the *.ctl file.
      --reference_protein_file REFERENCE_PROTEIN_FILE
                            If a given DNA sequence can be translated in more than one 
                            forward frame without stop codons, choose the translation 
                            that best matches the corresponding reference protein 
                            provided in this fasta file. The fasta headers must follow 
                            the naming convention: >TaxonID-geneName
      --maximum_distance FLOAT
                            When comparing candidate DNA translation frames to a 
                            reference protein, the maximum distance allowed between the 
                            translated frame and the reference sequence for any candidate 
                            translation frame to be selected. Useful to filter out 
                            sequences with frameshifts that do NOT introduce stop codons. 
                            0.0 means identical sequences, 1.0 means completely different 
                            sequences. Default is 0.5
      --filter_by_length_percentage FLOAT
                            If more than one representative sequence is present for a 
                            given gene, filter out sequences shorter than this percentage 
                            of the longest gene sequence length. Default is 0.0 (all 
                            sequences retained).
      --keep_low_complexity_sequences
                            Keep sequences that contain regions of low complexity, as 
                            identified by the command "hybpiper check_targetfile". Default 
                            is to remove these sequences.
      --alignments          Create per-gene alignments from the final fixed/filtered 
                            target file sequences. Note that DNA sequences will be 
                            translated prior to alignment.
      --concurrent_alignments INTEGER
                            Number of alignments to run concurrently. Default is 1.
      --threads_per_concurrent_alignment INTEGER
                            Number of threads to run each concurrent alignment with. 
                            Default is 1.
      --write_all_fasta_files
                            If provided, *.fasta files will be written for sequences 
                            removed from the fixed/filtered target file, according to 
                            filtering categories (length threshold, low-complexity 
                            regions, etc.). By default, these files will not be written.
      --verbose_logging     If supplied, enable verbose logging. NOTE: this will 
                            increase the size of the log files.
      --run_profiler        If supplied, run the subcommand using cProfile. Saves a 
                            *.csv file of results

  #######################################################################################
  """.stripIndent()
}

// Check for `--help` parameter if running `-entry fix_targetfile`:
if (workflow.commandLine.contains('-entry fix_targetfile') && params.remove('help') ) {
  fix_targetfile_help()
  exit 0
} 

// Check minimal input provided if `-entry fix_targetfile`:
if (workflow.commandLine.contains('-entry fix_targetfile')  && 
(!params.targetfile_dna && !params.targetfile_aa)) {
  println(
  """
  ERROR: Parameter "-entry fix_targetfile" provided, but no target file provided!
  Please provide your target file using the "--targetfile_dna" or "--targetfile_aa" parameter.
  """
  )
} 

if (workflow.commandLine.contains('-entry fix_targetfile')  && 
(!params.control_file)) {
  println(
  """
  ERROR: Parameter "-entry fix_targetfile" provided, but no control file provided!
  Please provide your control file using the "--control_file" parameter.
  """
  )
  exit 0
} 

// Create `hybpiper fix_targetfile` command string:
def fix_targetfile_command_list = []

fix_targetfile_command_list << "${params.control_file}"

if (params.targetfile_dna) {
  File target_file = new File(params.targetfile_dna)
  target_file_basename = target_file.getName()
  fix_targetfile_command_list << "--targetfile_dna ${target_file_basename}"
  }
if (params.targetfile_aa) {
  File target_file = new File(params.targetfile_aa)
  target_file_basename = target_file.getName()
  fix_targetfile_command_list << "--targetfile_aa ${target_file_basename}"
  }
if (params.no_terminal_stop_codons) {
  fix_targetfile_command_list << "--no_terminal_stop_codons"
  } 
if (params.allow_gene_removal) {
  fix_targetfile_command_list << "--allow_gene_removal"
  }
if (params.reference_protein_file) {
  fix_targetfile_command_list << "--reference_protein_file ${params.reference_protein_file}"
  }
if (params.maximum_distance) {
  fix_targetfile_command_list << "--maximum_distance ${params.maximum_distance}"
  }
if (params.filter_by_length_percentage) {
  fix_targetfile_command_list << "--filter_by_length_percentage ${params.filter_by_length_percentage}"
  }
if (params.keep_low_complexity_sequences) {
  fix_targetfile_command_list << "--keep_low_complexity_sequences ${params.keep_low_complexity_sequences}"
  }
if (params.alignments) {
  fix_targetfile_command_list << "--alignments"
  }
if (params.concurrent_alignments) {
  fix_targetfile_command_list << "--concurrent_alignments ${params.concurrent_alignments}"
  }
if (params.threads_per_concurrent_alignment) {
  fix_targetfile_command_list << "--threads_per_concurrent_alignment ${params.threads_per_concurrent_alignment}"
  }
if (params.write_all_fasta_files) {
  fix_targetfile_command_list << "--write_all_fasta_files"
  }
if (params.verbose_logging) {
  fix_targetfile_command_list << "--verbose_logging"
  }
if (params.run_profiler) {
  fix_targetfile_command_list << "--run_profiler ${params.run_profiler}"
  } 

// Workflow to run the fix_targetfile_main workflow with target_file_ch as input:
workflow fix_targetfile {
    fix_targetfile_main( target_file_ch, control_file_ch )
}

// Workflow to run the FIX_TARGETFILE process
workflow fix_targetfile_main {
    take: target_file
    take: control_file
    main:
        FIX_TARGETFILE(target_file, control_file)
        FIX_TARGETFILE.out.fix_results.view()
}


/////////////////////////////////////////////////////////////////////////////////////////
// HYBPIPER ASSEMBLE
/////////////////////////////////////////////////////////////////////////////////////////

def assemble_help() {
    log.info """

    Usage:
    The typical command for running the pipeline is as follows:

    nextflow run hybpiper.nf \
    -c hybpiper.config \
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
                                  standard_singularity, slurm_singularity, conda, 
                                  conda_slurm

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

      Options for step: map_reads:

      --bwa                       Use BWA to search reads for hits to target. Requires
                                  BWA and a target file that is nucleotides!

      --diamond                   Use DIAMOND instead of BLASTx

      --diamond_sensitivity       {mid-sensitive,sensitive,more-sensitive,very-sensitive,
                                   ultra-sensitive}
                                  Use the provided sensitivity for DIAMOND searches. 

      --evalue                    e-value threshold for blastx/DIAMOND hits, default: 0.0001

      --max_target_seqs           Max target seqs to save in BLASTx search, default is 10


      Options for step: distribute_reads:

      --distribute_low_mem        Distributing and writing reads to individual gene 
                                  directories  will be 40-50 percent slower, but can use 
                                  less memory/RAM with large input files

      
      Options for step: assemble_reads:

      --cov_cutoff <int>          Coverage cutoff to pass to the SPAdes assembler. 
                                  Default is 8

      --single_cell_assembly      Run SPAdes assemblies in MDA (single-cell) mode. 
                                  Default is False

      --kvals                     Values of k for SPAdes assemblies. SPAdes needs to be 
                                  compiled to handle larger k-values! Default is 
                                  auto-detection by SPAdes

      --merged                    Merge forward and reverse reads, and run SPAdes 
                                  assembly with merged and unmerged (the latter 
                                  in interleaved format) data. Default is off

      --timeout_assemble_reads    Kill long-running gene assemblies if they take longer 
                                  than X percent of average.


      Options for step: extract_contigs:

      --not_protein_coding        If provided, extract sequences from SPAdes contigs using 
                                  BLASTn rather than Exonerate (step: extract_contigs)

      --extract_contigs_blast_task
                                  {blastn,blastn-short,megablast,dc-megablast}
                                  Task to use for BLASTn searches during the extract_contigs 
                                  step of the assembly pipeline. See 
                                  https://www.ncbi.nlm.nih.gov/books/NBK569839/ for a 
                                  description of tasks. Default is: blastn

      --extract_contigs_blast_evalue
                                  Expectation value (E) threshold for saving hits. 
                                  Default is: 10
                                  
      --extract_contigs_blast_word_size
                                  Word size for wordfinder algorithm (length of best perfect 
                                  match)
      
      --extract_contigs_blast_gapopen
                                  Cost to open a gap

      --extract_contigs_blast_gapextend
                                  Cost to extend a gap

      --extract_contigs_blast_penalty
                                  Penalty for a nucleotide mismatch

      --extract_contigs_blast_reward
                                  Reward for a nucleotide match

      --extract_contigs_blast_perc_identity
                                  Percent identity. Can be used as a pre-filter at the BLASTn 
                                  stage, followed by --thresh (see below)

      --extract_contigs_blast_max_target_seqs
                                  Maximum number of aligned sequences to keep (value of 5 or 
                                  more is recommended). Default is: 500

      --thresh                    Percent identity threshold for retaining Exonerate/BLASTn
                                  hits. Default is 55, but increase this if you are worried 
                                  about contaminant sequences. Exonerate hit identity is 
                                  calculated using amino-acids, BLASTn hit identity is calculated 
                                  using nucleotides

      --paralog_min_length_percentage <decimal> 
                                  Minimum length percentage of a SPAdes contig vs 
                                  reference protein query for a paralog warning to be 
                                  generated and a putative paralog contig to be 
                                  recovered. Default is 0.75 

      --depth_multiplier          Assign a long paralog as the "main" sequence if it 
                                  has a coverage depth <depth_multiplier> times all 
                                  other long paralogs. Set to zero to not use depth. 
                                  Default is 10

      --target                    Use the target file sequence with this taxon name in 
                                  Exonerate searches for each gene. Other targets for 
                                  that gene will be used only for read sorting. Can be a 
                                  tab-delimited file (one <gene>\\t<taxon_name> per line) 
                                  or a single taxon name

      --exclude                   Do not use any sequence with the specified taxon name 
                                  string in Exonerate searches. Sequenced from this 
                                  taxon will still be used for read sorting

      --timeout_extract_contigs Kill long-running processes if they take longer than 
                                  X seconds. Default is 120

      --no_stitched_contig        Do not create stitched contigs; use longest Exonerate 
                                  hit only. Default is off

      --no_pad_stitched_contig_gaps_with_n
                                  When constructing stitched contigs, do not pad any gaps 
                                  between hits (with respect to the "best" protein reference) 
                                  with a number of Ns corresponding to the reference gap multiplied 
                                  by 3 (Exonerate) or reference gap (BLASTn). Default is: True.

      --chimeric_stitched_contig_check
                                  Attempt to determine whether a stitched contig is a potential 
                                  chimera of contigs from multiple paralogs. Default is: False

      --bbmap_memory <int>        MB memory (RAM) to use for bbmap.sh if a chimera check is 
                                  performed during step extract_contigs. Default: is 1000

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

      --trim_hit_sliding_window_size <int>
                                  Size of the sliding window (amino acids for Exonerate, 
                                  nucleotides for BLASTn) when trimming hit termini. 
                                  Default is: 5 (Exonerate) or 15 (BLASTn)

      --exonerate_hit_sliding_window_thresh <int>
                                  Percentage similarity threshold for the sliding window 
                                  (amino acids for Exonerate, nucleotides for BLASTn) 
                                  when trimming hit termini. Default is: 75 (Exonerate) 
                                  or 65 (BLASTn)

      --exonerate_skip_hits_with_frameshifts
                                  Skip Exonerate hits where the SPAdes sequence contains 
                                  a frameshift. See: 
                                  https://github.com/mossmatters/HybPiper/wiki/Troubleshooting-common-
                                  issues,-and-recommendations#42-hits-where-the-spades-contig-contains-frameshifts. 
                                  Default is: False

      --exonerate_skip_hits_with_internal_stop_codons
                                  Skip Exonerate hits where the SPAdes sequence contains an 
                                  internal in-frame stop codon. See:
                                  https://github.com/mossmatters/HybPiper/wiki/Troubleshooting,-common-
                                  issues,-and-recommendations#31-sequences-containing-stop-codons. 
                                  A single terminal stop codon is allowed, but see option 
                                  "--exonerate_skip_hits_with_terminal_stop_codons" below. Default is: False.

      --exonerate_skip_hits_with_terminal_stop_codons
                                  Skip Exonerate hits where the SPAdes sequence contains a single 
                                  terminal stop codon. Only applies when option 
                                  "--exonerate_skip_hits_with_internal_stop_codons" is also provided. 
                                  Only use this flag if your target file exclusively contains 
                                  protein-coding genes with no stop codons included, and you would like 
                                  to prevent any in-frame stop codons in the output sequences. 
                                  Default is: False.

      --exonerate_refine_full
                                  Run Exonerate searches using the parameter "--refine full". 
                                  Default is: False.
      
      --no_intronerate            Do not run intronerate to recover fasta files for supercontigs 
                                  with introns (if present), and introns-only. If this flag is used, 
                                  fasta files in `subfolders 09_sequences_intron` and 
                                  `10_sequences_supercontig` will be empty

      --no_padding_supercontigs   If Intronerate is run, and a supercontig is created 
                                  by concatenating multiple SPAdes contigs, do not add 
                                  10 "N" characters between contig joins. By default, 
                                  Ns will be added

      --keep_intermediate_files   Keep all intermediate files and logs, which can be 
                                  useful for debugging. Default action is to delete 
                                  them, which greatly reduces the total file number

      --verbose_logging           If supplied, enable verbose login. NOTE: this can 
                                  increase the size of the log files by an order of 
                                  magnitude

      --compress_sample_folder
                                  Tarball and compress the sample folder after assembly 
                                  has completed (<sample_name>.tar.gz). Default is: False

      #######################################################################################
      ####################### hybpiper paralog_retriever options: ###########################
      #######################################################################################

      --paralogs_list_threshold_percentage PARALOGS_LIST_THRESHOLD_PERCENTAGE
                                  Percent of total number of samples and genes that must
                                  have paralog warnings to be reported in the
                                  <genes_with_paralogs.txt> report file. The default is
                                  0.0, meaning that all genes and samples with at least
                                  one paralog warning will be reported
      --figure_length FIGURE_LENGTH
                                  Length dimension (in inches) for the output heatmap
                                  file. Default is automatically calculated based on the
                                  number of genes
      --figure_height FIGURE_HEIGHT
                                  Height dimension (in inches) for the output heatmap
                                  file. Default is automatically calculated based on the
                                  number of samples
      --sample_text_size SAMPLE_TEXT_SIZE
                                  Size (in points) for the sample text labels in the
                                  output heatmap file. Default is automatically
                                  calculated based on the number of samples
      --gene_text_size GENE_TEXT_SIZE
                                  Size (in points) for the gene text labels in the
                                  output heatmap file. Default is automatically
                                  calculated based on the number of genes
      --heatmap_filetype {png,pdf,eps,tiff,svg}
                                  File type to save the output heatmap image as. Default
                                  is png
      --heatmap_dpi HEATMAP_DPI
                                  Dots per inch (DPI) for the output heatmap image.
                                  Default is 300

      #######################################################################################
      ######################## hybpiper recovery_heatmap options: ###########################
      #######################################################################################

      --figure_length FIGURE_LENGTH
                                  Length dimension (in inches) for the output heatmap
                                  file. Default is automatically calculated based on the
                                  number of genes
      --figure_height FIGURE_HEIGHT
                                  Height dimension (in inches) for the output heatmap
                                  file. Default is automatically calculated based on the
                                  number of samples
      --sample_text_size SAMPLE_TEXT_SIZE
                                  Size (in points) for the sample text labels in the
                                  output heatmap file. Default is automatically
                                  calculated based on the number of samples
      --gene_text_size GENE_TEXT_SIZE
                                  Size (in points) for the gene text labels in the
                                  output heatmap file. Default is automatically
                                  calculated based on the number of genes
      --heatmap_filetype {png,pdf,eps,tiff,svg}
                                  File type to save the output heatmap image as. Default
                                  is *.png
      --heatmap_dpi HEATMAP_DPI
                                  Dot per inch (DPI) for the output heatmap image.
                                  Default is 150

    """.stripIndent()
}


// Check for `--help` parameter if running `-entry assemble`
if (workflow.commandLine.contains('-entry assemble') && 
params.help) {
  // if (params.help || !params.illumina_reads_directory || (!params.targetfile_dna && !params.targetfile_aa)) {
    assemble_help()
    exit 0
  }

// Check minimal input provided if `-entry assemble`, then check for incompatible parameters:
if (workflow.commandLine.contains('-entry assemble') && 
(!params.illumina_reads_directory || (!params.targetfile_dna && !params.targetfile_aa))) {

  println(
  """
  Please provide a folder of reads using the parameter "--illumina_reads_directory", and a \
target file using one of the parameters:

      --targetfile_dna 
      --targetfile_aa

  Use `--entry assemble --help` for more details.
  """
  )
  exit 0
  } 
  
if (params.paralog_min_length_percentage < 0 || params.paralog_min_length_percentage >1) {
    println(
      """
      The value for --paralog_min_length_percentage should be between 0 and 1. 
      Your value is ${params.paralog_min_length_percentage}
      """.stripIndent()
      )
    exit 0
  } 
  
if (params.single_end && params.paired_and_single) {
    println(
    """
    Please use --single_end OR --paired_and_single, not both!
    """
    )
    exit 0
  } 

if (params.paired_and_single && params.use_trimmomatic) {
    println(
    """
    Trimmomatic can't be used with paired plus single reads yet - 
    let me know if this would be useful!
    """.stripIndent()
    )
    exit 0
  } 
  
if (params.targetfile_dna && params.targetfile_aa) {
    println(
    """
    Please use --targetfile_dna OR --targetfile_aa, not both!
    """
    )
    exit 0
  } 
  
if (params.targetfile_aa && params.bwa) {
    println(
    """
    You can not use BWA with a target file containing protein sequences.
    Please use BLASTx or DIAMOND, or provide a target file with nucleotide sequences.
    """
    )
    exit 0
  }


// Create `hybpiper assemble` command string:
def command_list = []

if (params.targetfile_dna) {
  File target_file = new File(params.targetfile_dna)
  target_file_basename = target_file.getName()
  command_list << "--targetfile_dna ${target_file_basename}"
  }
if (params.targetfile_aa) {
  File target_file = new File(params.targetfile_aa)
  target_file_basename = target_file.getName()
  command_list << "--targetfile_aa ${target_file_basename}"
  }

// Options for step: map_reads:
if (params.bwa) {
  command_list << "--bwa"
  }
if (params.diamond) {
  command_list << "--diamond"
  }
if (params.diamond_sensitivity) {
  command_list << "--diamond_sensitivity ${params.diamond_sensitivity}"
  }
if (params.evalue) {
  command_list << "--evalue ${params.evalue}"
  }
if (params.max_target_seqs) {
  command_list << "--max_target_seqs ${params.max_target_seqs}"
  }

// Options for step: distribute_reads:
if (params.distribute_low_mem) {
  command_list << "--distribute_low_mem"
  }

// Options for step: assemble_reads:
if (params.cov_cutoff) {
  command_list << "--cov_cutoff ${params.cov_cutoff}"
  }
if (params.single_cell_assembly) {
  command_list << "--single_cell_assembly"
  }
if (params.kvals) {
  command_list << "--kvals ${params.kvals}"
  }
if (params.merged) {
  command_list << "--merged"
  }
if (params.timeout_assemble_reads) {
  command_list << "--timeout_assemble_reads ${params.timeout_assemble_reads}"
  }

command_list << "--no_spades_eta"  // Always use this parameter

// Options for step: extract_contigs:
if (params.not_protein_coding) {
  command_list << "--not_protein_coding"
  }
if (params.extract_contigs_blast_task) {
  command_list << "--extract_contigs_blast_task ${params.extract_contigs_blast_task}"
  }
if (params.extract_contigs_blast_evalue) {
  command_list << "--extract_contigs_blast_evalue ${params.extract_contigs_blast_evalue}"
  }
if (params.extract_contigs_blast_word_size) {
  command_list << "--extract_contigs_blast_word_size ${params.extract_contigs_blast_word_size}"
  }
if (params.extract_contigs_blast_gapopen) {
  command_list << "--extract_contigs_blast_gapopen ${params.extract_contigs_blast_gapopen}"
  }
if (params.extract_contigs_blast_gapextend) {
  command_list << "--extract_contigs_blast_gapextend ${params.extract_contigs_blast_gapextend}"
  }
if (params.extract_contigs_blast_penalty) {
  command_list << "--extract_contigs_blast_penalty ${params.extract_contigs_blast_penalty}"
  }
if (params.extract_contigs_blast_reward) {
  command_list << "--extract_contigs_blast_reward ${params.extract_contigs_blast_reward}"
  }
if (params.extract_contigs_blast_perc_identity) {
  command_list << "--extract_contigs_blast_perc_identity ${params.extract_contigs_blast_perc_identity}"
  }
if (params.extract_contigs_blast_max_target_seqs) {
  command_list << "--extract_contigs_blast_max_target_seqs ${params.extract_contigs_blast_max_target_seqs}"
  }
if (params.thresh) {
  command_list << "--thresh ${params.thresh}"
  }
if (params.paralog_min_length_percentage) {
  command_list << "--paralog_min_length_percentage ${params.paralog_min_length_percentage}"
  }
if (params.depth_multiplier) {
  command_list << "--depth_multiplier ${params.depth_multiplier}"
  }
if (params.target) {
  command_list << "--target ${params.target}"
  }
if (params.exclude) {
  command_list << "--exclude ${params.exclude}"
  } 
if (params.timeout_extract_contigs) {
  command_list << "--timeout_extract_contigs ${params.timeout_extract_contigs}"
  }
if (params.no_stitched_contig) {
  command_list << "--no_stitched_contig"
  }
if (params.no_pad_stitched_contig_gaps_with_n) {
  command_list << "--no_pad_stitched_contig_gaps_with_n"
  }
if (params.chimeric_stitched_contig_check) {
  command_list << "--chimeric_stitched_contig_check"
  }
if (params.bbmap_memory) {
  command_list << "--bbmap_memory ${params.bbmap_memory}"
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
if (params.trim_hit_sliding_window_size) {
  command_list << "--trim_hit_sliding_window_size ${params.trim_hit_sliding_window_size}"
  }
if (params.trim_hit_sliding_window_thresh) {
  command_list << "--trim_hit_sliding_window_thresh ${params.trim_hit_sliding_window_thresh}"
  }
if (params.exonerate_skip_hits_with_frameshifts) {
  command_list << "--exonerate_skip_hits_with_frameshifts ${params.exonerate_skip_hits_with_frameshifts}"
  }
if (params.exonerate_skip_hits_with_internal_stop_codons) {
  command_list << "--exonerate_skip_hits_with_internal_stop_codons ${params.exonerate_skip_hits_with_internal_stop_codons}"
  }
if (params.exonerate_skip_hits_with_terminal_stop_codons) {
  command_list << "--exonerate_skip_hits_with_terminal_stop_codons ${params.exonerate_skip_hits_with_terminal_stop_codons}"
  }
if (params.exonerate_refine_full) {
  command_list << "--exonerate_refine_full"
  }
if (params.no_intronerate) {
  command_list << "--no_intronerate"
  }
if (params.no_padding_supercontigs) {
  command_list << "--no_padding_supercontigs"
  }

// General `hybpiper assemble` pipeline options:
if (params.compress_sample_folder) {
  command_list << "--compress_sample_folder"
  }
if (params.skip_targetfile_checks) {
  command_list << "--skip_targetfile_checks"
  }
if (params.keep_intermediate_files) {
  command_list << "--keep_intermediate_files"
  }

if (params.verbose_logging) {
  command_list << "--verbose_logging"
  }
if (params.run_profiler) {
  command_list << "--run_profiler ${params.run_profiler}"
  }


// Create `hybpiper recovery_heatmap` command string:
def recovery_heatmap_command_list = []

if (params.figure_length) {
  recovery_heatmap_command_list << "--figure_length ${params.figure_length}"
  }
if (params.figure_height) {
  recovery_heatmap_command_list << "--figure_height ${params.figure_height}"
  }
if (params.sample_text_size) {
  recovery_heatmap_command_list << "--sample_text_size ${params.sample_text_size}"
  }
if (params.gene_text_size) {
  recovery_heatmap_command_list << "--gene_text_size ${params.gene_text_size}"
  }
if (params.heatmap_filetype) {
  recovery_heatmap_command_list << "--heatmap_filetype ${params.heatmap_filetype}"
  }
if (params.heatmap_dpi) {
  recovery_heatmap_command_list << "--heatmap_dpi ${params.heatmap_dpi}"
  }


// Create `hybpiper paralog_retriever` command string:
def paralog_retriever_command_list = []

if (params.paralogs_list_threshold_percentage) {
  paralog_retriever_command_list << "--paralogs_list_threshold_percentage ${params.paralogs_list_threshold_percentage}"
  }
if (params.figure_length) {
  paralog_retriever_command_list << "--figure_length ${params.figure_length}"
  }
if (params.figure_height) {
  paralog_retriever_command_list << "--figure_height ${params.figure_height}"
  }
if (params.sample_text_size) {
  paralog_retriever_command_list << "--sample_text_size ${params.sample_text_size}"
  }
if (params.gene_text_size) {
  paralog_retriever_command_list << "--gene_text_size ${params.gene_text_size}"
  }
if (params.heatmap_filetype) {
  paralog_retriever_command_list << "--heatmap_filetype ${params.heatmap_filetype}"
  }
if (params.heatmap_dpi) {
  paralog_retriever_command_list << "--heatmap_dpi ${params.heatmap_dpi}"
  }



// Workflow to run the assemble_main workflow with target_file_ch as input:
workflow assemble {
    assemble_main( target_file_ch )
}

// Workflow to run the `hybpiper assemble` pipeline:
workflow assemble_main {
    take: target_file
    main:

      //  Create 'namelist.txt' file and associated channel:
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


      //  Illumina reads channel:
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

      // Set up input channels for `hybpiper assemble`:
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

      // Run `hybpiper assemble`:
      ASSEMBLE_PAIRED_AND_SINGLE_END( 
        target_file_ch, 
        assemble_with_unpaired_input_ch 
        )
      ASSEMBLE_PAIRED_END( 
        target_file_ch, 
        assemble_no_unpaired_input_ch 
        )
      ASSEMBLE_SINGLE_END( 
        target_file_ch, 
        assemble_with_single_end_only_input_ch 
        )

      // Run `hybpiper stats`:
      SUMMARY_STATS( 
        ASSEMBLE_PAIRED_AND_SINGLE_END.out.assemble_with_unPaired_ch.collect()
        .mix(ASSEMBLE_PAIRED_AND_SINGLE_END.out.assemble_with_unPaired_gz_ch).collect()
        .mix(ASSEMBLE_PAIRED_END.out.assemble_ch).collect()
        .mix(ASSEMBLE_PAIRED_END.out.assemble_gz_ch).collect()
        .mix(ASSEMBLE_SINGLE_END.out.assemble_with_single_end_ch).collect() 
        .mix(ASSEMBLE_SINGLE_END.out.assemble_with_single_end_gz_ch).collect(), 
        target_file_ch, 
        namelist_ch 
        ) 

      // Run `hybpiper recovery_heatmap`: 
      VISUALISE( SUMMARY_STATS.out.seq_lengths_file ) 

      // Run hybpiper `retrieve_sequences` for all sequence types:
      RETRIEVE_SEQUENCES( 
        ASSEMBLE_PAIRED_AND_SINGLE_END.out.assemble_with_unPaired_ch.collect()
        .mix(ASSEMBLE_PAIRED_AND_SINGLE_END.out.assemble_with_unPaired_gz_ch).collect()
        .mix(ASSEMBLE_PAIRED_END.out.assemble_ch).collect()
        .mix(ASSEMBLE_PAIRED_END.out.assemble_gz_ch).collect()
        .mix(ASSEMBLE_SINGLE_END.out.assemble_with_single_end_ch).collect() 
        .mix(ASSEMBLE_SINGLE_END.out.assemble_with_single_end_gz_ch).collect(), 
        target_file_ch, 
        namelist_ch 
        )

      // Run `hybpiper paralog_retriever`: 
      PARALOG_RETRIEVER( 
        ASSEMBLE_PAIRED_AND_SINGLE_END.out.assemble_with_unPaired_ch.collect()
        .mix(ASSEMBLE_PAIRED_AND_SINGLE_END.out.assemble_with_unPaired_gz_ch).collect()
        .mix(ASSEMBLE_PAIRED_END.out.assemble_ch).collect()
        .mix(ASSEMBLE_PAIRED_END.out.assemble_gz_ch).collect()
        .mix(ASSEMBLE_SINGLE_END.out.assemble_with_single_end_ch).collect() 
        .mix(ASSEMBLE_SINGLE_END.out.assemble_with_single_end_gz_ch).collect(), 
        namelist_ch, 
        target_file_ch )
      
}


/////////////////////////////////////////////////////////////////////////////////////////
//  DEFINE DSL2 PROCESSES
/////////////////////////////////////////////////////////////////////////////////////////


process CHECK_TARGETFILE {
  /*
  Run `hybpiper check_targetfile` command.
  */

  publishDir "${workflow.launchDir}", mode: 'copy', pattern: "check_targetfile_report.txt"
  publishDir "${workflow.launchDir}", mode: 'copy', pattern: "fix_targetfile_*.ctl"

  input:
    path(target_file)

  output:
    stdout emit: check_results
    path("check_targetfile_report.txt")
    path("fix_targetfile_*.ctl")

  script:
  check_targetfile_assemble_command = "hybpiper check_targetfile " + check_targetfile_command_list.join(' ') + "| tee check_targetfile_report.txt"

    """
    echo "Executing command: ${check_targetfile_assemble_command}"
    ${check_targetfile_assemble_command}
    """
}


process FIX_TARGETFILE {
  /*
  Run `hybpiper fix_targetfile` command.
  */

  publishDir "${workflow.launchDir}", mode: 'copy', pattern: "fix_targetfile_report.txt"
  publishDir "${workflow.launchDir}", mode: 'copy', pattern: "fix_targetfile_report.tsv"
  publishDir "${workflow.launchDir}", mode: 'copy', pattern: "fix_targetfile_*.log"
  publishDir "${workflow.launchDir}", mode: 'copy', pattern: "fix_targetfile_alignments"
  publishDir "${workflow.launchDir}", mode: 'copy', pattern: "fix_targetfile_additional_sequence_files"
  publishDir "${workflow.launchDir}", mode: 'copy', pattern: "*fixed*"
  

  input:
    path(target_file)
    path(control_file)

  output:
    stdout emit: fix_results
    path("fix_targetfile_report.txt")
    path("fix_targetfile_report.tsv")
    path("*fixed*")
    path("fix_targetfile_alignments"), optional: true
    path("fix_targetfile_additional_sequence_files"), optional: true

  script:
  fix_targetfile_assemble_command = "hybpiper fix_targetfile " + fix_targetfile_command_list.join(' ') + " | tee fix_targetfile_report.txt"

    """
    echo "Executing command: ${fix_targetfile_assemble_command}"
    ${fix_targetfile_assemble_command}
    """
}


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
  publishDir "${params.outdir}/04_processed_sample_directories", mode: 'copy', pattern: "${pair_id}"
  publishDir "${params.outdir}/04_processed_sample_directories", mode: 'copy', pattern: "${pair_id}.tar.gz"
  publishDir "${params.outdir}/06_summary_stats", mode: 'copy', pattern: "${pair_id}/${pair_id}_genes_with_stitched_contig.csv"
  publishDir "${params.outdir}/06_summary_stats", mode: 'copy', pattern: "${pair_id}/${pair_id}_genes_with_long_paralog_warnings.txt"
  publishDir "${params.outdir}/06_summary_stats", mode: 'copy', pattern: "${pair_id}/${pair_id}_genes_with_paralog_warnings_by_contig_depth.csv"

  if (params.num_forks) {
    maxForks params.num_forks
  }

  when:
    params.single_end

  input:
    path(target_file)
    tuple val(prefix), path(reads_single)

  output:
    path("${prefix}"), emit: assemble_with_single_end_ch, optional: true
    path("${prefix}.tar.gz"), emit: assemble_with_single_end_gz_ch, optional: true
    path("${pair_id}/${pair_id}_genes_with_stitched_contig.csv"), optional: true
    path("${pair_id}/${pair_id}_genes_with_long_paralog_warnings.txt"), optional: true
    path("${pair_id}/${pair_id}_genes_with_paralog_warnings_by_contig_depth.csv"), optional: true

  script:
    assemble_command = "hybpiper assemble -r ${reads_single} --prefix ${prefix} --cpu ${task.cpus} " + command_list.join(' ')

    """
    echo ${assemble_command}
    ${assemble_command}
    """
}


process ASSEMBLE_PAIRED_AND_SINGLE_END {
  /*
  Run the `hybpiper assemble` command for input files: [R1, R1, R1-R2_unpaired]
  */

  //echo true
  label 'in_container'
  publishDir "${params.outdir}/04_processed_sample_directories", mode: 'copy', pattern: "${pair_id}"
  publishDir "${params.outdir}/04_processed_sample_directories", mode: 'copy', pattern: "${pair_id}.tar.gz"
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
    path("${pair_id}"), emit: assemble_with_unPaired_ch, optional: true
    path("${pair_id}.tar.gz"), emit: assemble_with_unPaired_gz_ch, optional: true
    path("${pair_id}/${pair_id}_genes_with_stitched_contig.csv"), optional: true
    path("${pair_id}/${pair_id}_genes_derived_from_putative_chimeric_stitched_contig.csv"), optional: true
    path("${pair_id}/${pair_id}_genes_with_long_paralog_warnings.txt"), optional: true
    path("${pair_id}/${pair_id}_genes_with_paralog_warnings_by_contig_depth.csv"), optional: true

  script:
    assemble_command = "hybpiper assemble -r ${reads_R1} ${reads_R2} --unpaired ${reads_unpaired} --prefix ${pair_id} --cpu ${task.cpus} " + command_list.join(' ')

    """
    echo ${assemble_command}
    ${assemble_command}
    """
  } 


process ASSEMBLE_PAIRED_END {
  /*
  Run the `hybpiper assemble` command for input files: [R1, R1]
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}/04_processed_sample_directories", mode: 'copy', pattern: "${pair_id}"
  publishDir "${params.outdir}/04_processed_sample_directories", mode: 'copy', pattern: "${pair_id}.tar.gz"
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
    path("${pair_id}"), emit: assemble_ch, optional: true
    path("${pair_id}.tar.gz"), emit: assemble_gz_ch, optional: true
    path("${pair_id}/${pair_id}_genes_with_stitched_contig.csv"), optional: true
    path("${pair_id}/${pair_id}_genes_derived_from_putative_chimeric_stitched_contig.csv"), optional: true
    path("${pair_id}/${pair_id}_genes_with_long_paralog_warnings.txt"), optional: true
    path("${pair_id}/${pair_id}_genes_with_paralog_warnings_by_contig_depth.csv"), optional: true


  script:
    assemble_command = "hybpiper assemble -r ${reads_R1} ${reads_R2} --prefix ${pair_id} --cpu ${task.cpus} " + command_list.join(' ')

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
  File target_file = new File(params.targetfile_dna)
  target_file_basename = target_file.getName()
  """
  hybpiper stats --cpu ${task.cpus} -t_dna ${target_file_basename} gene ${namelist}
  """
  } else if (params.targetfile_aa) {
  File target_file = new File(params.targetfile_aa)
  target_file_basename = target_file.getName()
  """
  hybpiper stats --cpu ${task.cpus} -t_aa ${target_file_basename} gene ${namelist}
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
    heatmap_command = "hybpiper recovery_heatmap ${seq_lengths_file} " + recovery_heatmap_command_list.join(' ')

    """
    echo "Executing command: ${heatmap_command}"
    ${heatmap_command}
    """
}


process RETRIEVE_SEQUENCES {
  /*
  Run the hybpiper retrieve_sequences command for all sequence types.
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}/07_sequences_dna", mode: 'copy', pattern: "00_dna_seqs/*.FNA", saveAs: { filename -> file(filename).getName() }
  publishDir "${params.outdir}/08_sequences_aa", mode: 'copy', pattern: "01_aa_seqs/*.FAA", saveAs: { filename -> file(filename).getName() }
  publishDir "${params.outdir}/09_sequences_intron", mode: 'copy', pattern: "02_intron_seqs/*.fasta", saveAs: { filename -> file(filename).getName() }
  publishDir "${params.outdir}/10_sequences_supercontig", mode: 'copy', pattern: "03_supercontig_seqs/*.fasta", saveAs: { filename -> file(filename).getName() }

  input:
    path(assemble)
    path(target_file)
    path(namelist)

  output:
    path("00_dna_seqs/*.FNA")
    path("01_aa_seqs/*.FAA")
    path("02_intron_seqs/*.fasta")
    path("03_supercontig_seqs/*.fasta")

  script:
  if (params.targetfile_dna) {
    File target_file = new File(params.targetfile_dna)
    target_file_basename = target_file.getName()
    """
    hybpiper retrieve_sequences --cpu ${task.cpus} -t_dna ${target_file_basename} --sample_names ${namelist} --fasta_dir 00_dna_seqs dna
    hybpiper retrieve_sequences --cpu ${task.cpus} -t_dna ${target_file_basename} --sample_names ${namelist} --fasta_dir 01_aa_seqs aa
    hybpiper retrieve_sequences --cpu ${task.cpus} -t_dna ${target_file_basename} --sample_names ${namelist} --fasta_dir 02_intron_seqs intron
    hybpiper retrieve_sequences --cpu ${task.cpus} -t_dna ${target_file_basename} --sample_names ${namelist} --fasta_dir 03_supercontig_seqs supercontig
    """
  } else if (params.targetfile_aa) {
    File target_file = new File(params.targetfile_aa)
    target_file_basename = target_file.getName()
    """
    hybpiper retrieve_sequences --cpu ${task.cpus} -t_aa ${target_file_basename} --sample_names ${namelist} --fasta_dir 00_dna_seqs dna
    hybpiper retrieve_sequences --cpu ${task.cpus} -t_aa ${target_file_basename} --sample_names ${namelist} --fasta_dir 01_aa_seqs aa
    hybpiper retrieve_sequences --cpu ${task.cpus} -t_aa ${target_file_basename} --sample_names ${namelist} --fasta_dir 02_intron_seqs intron
    hybpiper retrieve_sequences --cpu ${task.cpus} -t_aa ${target_file_basename} --sample_names ${namelist} --fasta_dir 03_supercontig_seqs supercontig
    """
  }
}


process PARALOG_RETRIEVER {
  /*
  Run hybpiper `paralog_retriever` command.
  */

  //echo true
  label 'in_container'
  publishDir "${params.outdir}/05_visualise", mode: 'copy', pattern: "paralog_heatmap.*"
  publishDir "${params.outdir}/11_paralogs", mode: 'copy', pattern: "paralogs_all/*paralogs_all.fasta", saveAs: { filename -> file(filename).getName() }
  publishDir "${params.outdir}/12_paralogs_no_chimeras", mode: 'copy', pattern: "paralogs_no_chimeras/*paralogs_no_chimeras.fasta", saveAs: { filename -> file(filename).getName() }
  publishDir "${params.outdir}/11_paralogs/logs", mode: 'copy', pattern: "*report*"

  input:
    path(assemble)
    path(namelist)
    path(target_file)


  output:
    path("paralogs_all/*paralogs_all.fasta")
    path("paralogs_no_chimeras/*paralogs_no_chimeras.fasta"), optional: true
    path("paralog_report.tsv")
    path("paralogs_above_threshold_report.txt") 
    path("paralog_heatmap.*")

  script:
    if (params.targetfile_dna) {
      File target_file = new File(params.targetfile_dna)
      target_file_basename = target_file.getName()
      paralog_command = "hybpiper paralog_retriever ${namelist} -t_dna ${target_file_basename} " + paralog_retriever_command_list.join(' ')

    } else if (params.targetfile_aa) {
      File target_file = new File(params.targetfile_aa)
      target_file_basename = target_file.getName()
      paralog_command = "hybpiper paralog_retriever ${namelist} --cpu ${task.cpus} -t_aa ${target_file_basename} " + paralog_retriever_command_list.join(' ')

      """
      echo "Executing command: ${paralog_command}"
      ${paralog_command}
      """
    }
}


/////////////////////////////////////////////////////////////////////////////////////////
// END OF SCRIPT  
/////////////////////////////////////////////////////////////////////////////////////////
