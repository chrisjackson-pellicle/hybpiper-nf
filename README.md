# HybPiper-RBGV

## Original HybPiper github, wiki and tutorial:

https://github.com/mossmatters/HybPiper/

https://github.com/mossmatters/HybPiper/wiki

https://github.com/mossmatters/HybPiper/wiki/Tutorial

## RBGB-modified HybPiper: containerised and pipelined using Singularity and Nextflow

To simplify running HybPiper, I’ve made a Singularity container containing the Linux distribution Ubuntu 18.04, containing all the scripts required to run the HybPiper pipeline (including modifications and bug fixes, see below for details), as well as all the dependencies (BioPython, BWA, BBmap [new requirement compared to default HybPiper], Exonerate, SPAdes, Samtools). The container is called `hybpiper_only.sif`.

To run the pipeline, I’ve made a Nextflow script that uses the software in the Singularity container. This script runs all HybPiper steps with a single command. The script is currently called `hybpiper_pipeline_v1_7.nf`. It comes with an associated config file called `hybpiper_v1.7.config`. The only input required is a folder of sequencing reads for your samples, and a target file in fasta format. The Nextflow pipeline will automatically generate the `namelist.txt` file, and will run all HybPiper scripts on each sample in parallel. The number of parallel processes running at any time, as well as computing resources given to each process (e.g. number of CPUs, amount of RAM etc) can be configured by the user by modifying a config file that comes with the Nextflow pipeline script. The pipeline can be run directly on your local computer, and on an HPC system submitting jobs via a scheduler (e.g. SLURM, PBS, etc).

## Name formatting of input read files
You will need to provide the `hybpiper_pipeline_v1_7.nf` script with a directory of forwards and reverse reads (and, optionally, a file of single reads) for each sample. See below for run command details. For the read files to be recognised by the script, they should be named according to the default convention:

    *_R1.fastq 
    *_R2.fastq
    *_single.fastq (optional - will be used if running with the flag `--unpaired`)

OR

    *_R1.fq 
    *_R2.fq
    *_single.fq (optional - will be used if running with the flag `--unpaired`)


It’s fine if there’s text after the R1/R1 and before the .fastq/.fq.
It’s fine if the input files are gzipped (i.e. suffix .gz). 

You can also specify the pattern used for file matching via the parameters `--read_pairs_pattern <pattern>` or  `--single_pattern <pattern>`. You can also just provide a folder of single-end reads (use the flag `--single_only` if you do)


## Running on Linux

Please see the Wiki entry [Running on Linux]

## Running on a Mac (macOS)

Please see the Wiki entry [Running on a Mac]

## Nextflow pipeline options and parameters

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

     --single_only                             Use when providing providing only a folder of single-end reads

     --outdir <directory_name>                 Specify the name of the pipeline results directory. Default is 'results'

     --read_pairs_pattern <pattern>            Provide a comma-separated read pair pattern for matching forwards and 
                                               reverse paired-end read-files. Default is 'R1,R2'

     --single_pattern <pattern>                Provide a pattern for matching single-end read files. Default is 'single'

     --use_blastx                              Use a protein target file and map reads to targets with BLASTx. Default 
                                               is a nucleotide target file and mapping of reads to targets using BWA

     --num_forks <int>                         Specify the number of parallel processes (e.g. concurrent runs of 
                                               reads.first.py) to run at any one time. Can be used to prevent Nextflow 
                                               from using all the threads/cpus on your machine. Default is to use the 
                                               maximum number possible

     --cov_cutoff <int>                        Coverage cutoff to pass to the SPAdes assembler. Default is 8

     --blastx_evalue <value>                   Evalue to pass to blastx when using blastx mapping (i.e. when the 
                                               --use_blastx flag is specified). Default is 1e-4

     --paralog_warning_min_len_percent <decimal>
                                               Minimum length percentage of a contig vs reference protein length for 
                                               a paralog warning to be generated and a putative paralog contig to be recovered.
                                               Default is 0.75

     --translate_target_file_for_blastx        Translate a nucleotide target file. If set, the --use_blastx is set by default. Default is off

     --use_trimmomatic                         Trim forwards and reverse reads using Trimmomatic. Default is off

     --trimmomatic_leading_quality <int>       Cut bases off the start of a read, if below this threshold quality. Default is 3

     --trimmomatic_trailing_quality <int>      Cut bases off the end of a read, if below this threshold quality. Default 3

     --trimmomatic_min_length <int>            Drop a read if it is below this specified length. Default is 36

     --trimmomatic_sliding_window_size <int>   Size of the sliding window used by Trimmomatic; specifies the number of 
                                               bases to average across. Default is 4

     --trimmomatic_sliding_window_quality <int>
                                               Specifies the average quality required within the sliding window. Default is 20

     --run_intronerate                         Run intronerate.py to recover (hopefully) intron and supercontig sequences.
                                               Default is off, and so results subfolders 09_sequences_intron and 
                                               10_sequences_supercontig will be empty
                                               
                                               
## General notes

## Bug fixed and changes (WIP)

- Read mapping to detect putative chimeras
- Exonerate -refine full with fallback
- SPAdes SC mode? 
- Reporting when supercontigs are made for a sample/locus (i.e. when multiple SPAdes contigs are concatenated together) 
- Reporting when supercontig creation requires trimming of overlaps between concatenated contigs (cause by overlaps in Exonerate hits, likely due to the contigs originating from different paralogs rather than being exons of the same gene)
- Misc fixes
   - Fixed a ‘ZeroDivisionError: float division by zero’ error in hybpiper_stats.py, cause when no reads map to the target file (see function enrich_efficiency_bwa())
   - Range tests fix so that good contigs with same exonerate hits aren’t both thrown out (see function tuple_subsume())
   - 
- No supercontigs option i.e. just return the longest single SPAdes contig
- Merged option - optionally merge paired end reads and run SPAdes assembly with merged and unmerged - gives better results when your Illumina library fragment size is small enough for R1 and R2 reads to overlap
- .gz file input option - allows input read files to be compressed in .gz format
- Additional stats reporting

## Issues still to deal with

- SPAdes can fail during assembly of some kmer lengths in a manner that isn’t detected by HybPiper (looks like this is a bug in SPAdes), meaning that real contigs are missed for some loci.
- SPAdes can hang indefinitely at the kmer coverage model stage (see https://github.com/ablab/spades/issues/653), meaning the pipeline never finishes. Running SPAdes in single-cell mode (--SC) seems to fix this, but it can also create many more spurious contigs than running without single-cell mode, which can in turn sometimes cause errors in the loci sequence output.
- Detection of paralogs - it’s currently very rough, and will miss many paralogs.
- Intronerate.py can error out - needs further investigation but looks like a bug in the way Exonerate output is parsed, which throws its results into question even when it does run successfully. UPDATE v1.7: I think I’ve fixed the bug causing the error, but I still wouldn’t use the output of intronerate unless you’ve manually vetted AND QC’d all the sequences.
- How does Intronerate.py work with the changes I’ve made to supercontig assembly (i.e. trimming overlaps between contigs before they’re concatenated, to avoid the introduction of repeats in the output loci sequence)? UPDATE v1.7: I’ve now checked on this, and it shouldn’t be negatively affected. The Intronerate.py script uses the protein .FAA sequence output by the reads_first.py script as a query, and re-runs Exonerate using the full SPAdes contigs as subjects (those that returned Exonerate hits in the reads_first.py run, anyway). So, the changes I’ve made should just (theoretically) improve the quality of input protein for the Intronerate.py run. 
- If a target file sequence for a particular gene has a low-complexity area, many reads can map to it and it will be selected and translated as the protein reference for a particular sample/gene, and then used in Exonerate searches of SPAdes contigs. This protein reference may not be the most similar to the particular sample/gene, meaning that Exonerate fails to identify (or simply truncates) real contigs for this locus. 
- Paralog tests do not remove flagged contigs from downstream assembly, so they can still get stitched into a supercontig to form a likely chimera.






