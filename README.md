# HybPiper-RBGV

## Original HybPiper github, wiki and tutorial:

For an explanation of the general purpose of HybPiper, and the approach it takes to generate sequences for target loci, please see the excellent documentation, wiki and tutorial at the `mossmatters` repo:

- https://github.com/mossmatters/HybPiper/
- https://github.com/mossmatters/HybPiper/wiki
- https://github.com/mossmatters/HybPiper/wiki/Tutorial

## RBGB-modified HybPiper: containerised and pipelined using Singularity and Nextflow

To simplify running HybPiper, I’ve made a Singularity container containing the Linux distribution Ubuntu 18.04, containing all the scripts required to run the HybPiper pipeline (including modifications and bug fixes, see below for details), as well as all the dependencies (BioPython, BWA, BBmap [new requirement compared to default HybPiper], Exonerate, SPAdes, Samtools). The container is called `hybpiper_only.sif`.

To run the pipeline, I’ve made a Nextflow script that uses the software in the Singularity container. This script runs all HybPiper steps with a single command. The script is currently called `hybpiper_pipeline_v1_7.nf`. It comes with an associated config file called `hybpiper_v1.7.config`. The only input required is a folder of sequencing reads for your samples, and a target file in fasta format. The Nextflow pipeline will automatically generate the `namelist.txt` file, and will run all HybPiper scripts on each sample in parallel. The number of parallel processes running at any time, as well as computing resources given to each process (e.g. number of CPUs, amount of RAM etc) can be configured by the user by modifying a config file that comes with the Nextflow pipeline script. The pipeline can be run directly on your local computer, and on an HPC system submitting jobs via a scheduler (e.g. SLURM, PBS, etc).

## Name formatting of input read files
You will need to provide the `hybpiper_pipeline_v1_7.nf` pipeline with either:

a) A directory of forwards and reverse reads (and, optionally, a file of single reads) for each sample; or
b) A directory of single-end reads.

For the read files to be recognised by the script, they should be named according to the default convention:

    *_R1.fastq 
    *_R2.fastq
    *_single.fastq (optional - will be used if running with the flag `--unpaired`)

OR

    *_R1.fq 
    *_R2.fq
    *_single.fq (optional - will be used if running with the flag `--unpaired`)

**NOTE:**
- It’s fine if there’s text after the R1/R1 and before the .fastq/.fq.
- It’s fine if the input files are gzipped (i.e. suffix .gz). 
- You can also specify the pattern used for file matching via the parameters `--read_pairs_pattern <pattern>` or  `--single_pattern <pattern>`. You can also just provide a folder of single-end reads (use the flag `--single_only` if you do)


- The pipeline will concatenate samples that have been run on different lanes, e.g. read files:

      79678_LibID81729_HF7CKAFX2_TGAATGCC-TGTCTAGT_L001_R1.fastq
      79678_LibID81729_HF7CKAFX2_TGAATGCC-TGTCTAGT_L001_R2.fastq
      79678_LibID81729_HF7CKAFX2_TGAATGCC-TGTCTAGT_L002_R1.fastq
      79678_LibID81729_HF7CKAFX2_TGAATGCC-TGTCTAGT_L002_R2.fastq
      79679_LibID81730_HF7CKAFX2_GCAACTAT-TCGTTGAA_L001_R1.fastq
      79679_LibID81730_HF7CKAFX2_GCAACTAT-TCGTTGAA_L001_R2.fastq
      79679_LibID81730_HF7CKAFX2_GCAACTAT-TCGTTGAA_L002_R1.fastq
      79679_LibID81730_HF7CKAFX2_GCAACTAT-TCGTTGAA_L002_R2.fastq

  It does this by grouping forwards and reverse reads (or single reads if you’re providing a folder of single reads and using the `--single_only` flag) via the common prefix preceding the first underscore (‘_’). So, `79678` and `79679` in this case.


## Running on Linux

Please see the Wiki entry [Running on Linux][2]

## Running on a Mac (macOS)

Please see the Wiki entry [Running on a Mac][1]

## Nextflow pipeline options and parameters

Example run command:

    nextflow run hybpiper_pipeline_v1_7.nf -c hybpiper_v1_7.config --illumina_reads_directory reads_for_hybpiper --target_file Angiosperms353_targetSequences.fasta


Options:


    Mandatory arguments:

    ##############################################################################

    --illumina_reads_directory <directory>    Path to folder containing illumina read file(s)
  
    --target_file <file>                      File containing fasta sequences of target genes

    ##############################################################################

    Optional arguments:

     -profile <profile>                       Configuration profile to use. Can use multiple (comma separated)
                                              Available: standard (default), slurm

     --cleanup                                Run 'cleanup.py' for each gene directory after 'reads_first.py'

     --nosupercontigs                         Do not create supercontigs. Use longest Exonerate hit only. Default is off.

     --memory <int>                           Memory (RAM) amount in GB to use for bbmap.sh with exonerate_hits.py. 
                                              Default is 1 GB

     --discordant_reads_edit_distance <int>   Minimum number of base differences between one read of a read pair 
                                              vs the supercontig reference for a read pair to be flagged as discordant. 
                                              Default is 5

     --discordant_reads_cutoff <int>          Minimum number of discordant reads pairs required to flag a supercontigs 
                                              as a potential hybrid of contigs from multiple paralogs. Default is 5

     --merged                                 Merge forward and reverse reads, and run SPAdes assembly with merged and 
                                              unmerged (the latter in interleaved format) data. Default is off

     --paired_and_single                      Use when providing both paired R1 and R2 read files as well as a file of 
                                              single-end reads for each sample

     --single_only                            Use when providing providing only a folder of single-end reads

     --outdir <directory_name>                Specify the name of the pipeline results directory. Default is 'results'

     --read_pairs_pattern <pattern>           Provide a comma-separated read pair pattern for matching forwards and 
                                              reverse paired-end read-files. Default is 'R1,R2'

     --single_pattern <pattern>               Provide a pattern for matching single-end read files. Default is 'single'

     --use_blastx                             Use a protein target file and map reads to targets with BLASTx. Default 
                                              is a nucleotide target file and mapping of reads to targets using BWA

     --num_forks <int>                        Specify the number of parallel processes (e.g. concurrent runs of 
                                              reads.first.py) to run at any one time. Can be used to prevent Nextflow 
                                              from using all the threads/cpus on your machine. Default is to use the 
                                              maximum number possible

     --cov_cutoff <int>                       Coverage cutoff to pass to the SPAdes assembler. Default is 8

     --blastx_evalue <value>                  Evalue to pass to blastx when using blastx mapping (i.e. when the 
                                              --use_blastx flag is specified). Default is 1e-4

     --paralog_warning_min_len_percent <decimal>
                                              Minimum length percentage of a contig vs reference protein length for 
                                              a paralog warning to be generated and a putative paralog contig to be recovered.
                                              Default is 0.75

     --translate_target_file_for_blastx       Translate a nucleotide target file. If set, the --use_blastx is set by default. Default is off

     --use_trimmomatic                        Trim forwards and reverse reads using Trimmomatic. Default is off

     --trimmomatic_leading_quality <int>      Cut bases off the start of a read, if below this threshold quality. Default is 3

     --trimmomatic_trailing_quality <int>     Cut bases off the end of a read, if below this threshold quality. Default 3

     --trimmomatic_min_length <int>           Drop a read if it is below this specified length. Default is 36

     --trimmomatic_sliding_window_size <int>  Size of the sliding window used by Trimmomatic; specifies the number of 
                                              bases to average across. Default is 4

     --trimmomatic_sliding_window_quality <int>
                                              Specifies the average quality required within the sliding window. Default is 20

     --run_intronerate                        Run intronerate.py to recover (hopefully) intron and supercontig sequences.
                                              Default is off, and so results subfolders 09_sequences_intron and 
                                              10_sequences_supercontig will be empty
                                            
Please see the Wiki entry [Additional pipeline features and details][5] for further explanation of the parameters above, and general pipeline functionality.
             

## Output folders and files

If you're just after the unaligned `.fasta` files for each of your target genes, ready for concatenation or coalescent-based analyses, the two main output folders of interest are probably:

- `07_sequences_dna`
- `08_sequences_aa`

For a full explanation of output folders and files, please see the Wiki entry [Output folders and files][6]
             
## General notes

- Link to Chimera detection approach
- Description of overlapping Exonerate hits issue

## Bug fixes and changes (WIP)

Please see the Wiki entry [Bug fixed and changes][2]



## Issues still to deal with (WIP)

Please see the Wiki entry [Issues][4]




[1]:https://github.com/chrisjackson-pellicle/HybPiper-RBGV/wiki/Running-on-a-Mac-(macOS)-with-Vagrant "Link to Running on a Mac Wiki"
[2]:https://github.com/chrisjackson-pellicle/HybPiper-RBGV/wiki/Bug-fixes-and-changes-(WIP) "Link to bug fixes and changes Wiki"
[3]:https://github.com/chrisjackson-pellicle/HybPiper-RBGV/wiki/Running-on-Linux "Link to Running on Linux Wiki"
[4]:https://github.com/chrisjackson-pellicle/HybPiper-RBGV/wiki/Issues "Link to Issues Wiki"
[5]:https://github.com/chrisjackson-pellicle/HybPiper-RBGV/wiki/Additional-pipeline-features-and-details "Link to Additional pipeline features and details Wiki"
[6]:https://github.com/chrisjackson-pellicle/HybPiper-RBGV/wiki/Output-folders-and-files "Link to Output folders and files Wiki"
