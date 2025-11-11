# hybpiper-nf

## Original HybPiper github, wiki and tutorial:

For an explanation of the general purpose of [HybPiper][17], and the approach it takes to generate target locus sequences from sequence capture data, please see the documentation, wiki and tutorial at the `mossmatters` repo:

- https://github.com/mossmatters/HybPiper/
- https://github.com/mossmatters/HybPiper/wiki
- https://github.com/mossmatters/HybPiper/wiki/Tutorial

## hybpiper-nf: containerised and pipelined using Singularity and Nextflow

To simplify running HybPiper, I’ve provided a [Singularity][15] container based on the Linux distribution Ubuntu 22.04, with all the software required to run the HybPiper pipeline (including some additional functionality, see below for details), as well as all the dependencies ([BioPython][9], [BLAST][8], [BWA][11], [BBmap][12], [Exonerate][13], [SPAdes][10], [Samtools][14]). The container is called `hybpiper-paragone.sif`.

To run HybPiper using this container, I’ve provided a [Nextflow][16] pipeline that uses the software in the Singularity container. This pipeline runs all HybPiper steps with a single command. The pipeline script is called `hybpiper.nf`. It comes with an associated config file called `hybpiper.config`. The only input required is a folder of sequencing reads for your samples, and a target file in `.fasta` format. The Nextflow pipeline will automatically generate the `namelist.txt` file required by some of the HybPiper scripts, and will run all HybPiper scripts on each sample in parallel. It also includes an optional read-trimmming step to QC your reads prior to running HybPiper, using the software [Trimmomatic][7]. The number of parallel processes running at any time, as well as computing resources given to each process (e.g. number of CPUs, amount of RAM etc) can be configured by the user by modifying the provided config file. The pipeline can be run directly on your local computer, and on an HPC system submitting jobs via a scheduler (e.g. [SLURM][21], PBS, etc). 

## Input data

### Target file

Your **`target file`**, which can contain either nucleotide OR protein sequences (see below for corresponding pipeline options), should follow the formatting described for HybPiper [here][19]. Briefly, your target genes should be grouped and differentiated by a suffix in the fasta header, consisting of a dash followed by an ID unique to each gene, e.g.:

    >AJFN-4471
    AATGTTATACAGGATGAAGAGAAACTGAATACTGCAAACTCCGATTGGATGCGGAAATACAAAGGCT...
    >Ambtr-4471
    AGTGTTATTCAAGATGAAGATGTATTGTCGACAGCCAATGTGGATTGGATGCGGAAATATAAGGGCA...
    >Ambtr-4527
    GAGGAGCGGGTGATTGCCTTGGTCGTTGGTGGTGGGGGTAGAGAACATGCTCTATGCTATGCTTTGC...
    >Arath-4691
    GAGCTTGGATCTATCGCTTGCGCAGCTCTCTGTGCTTGCACTCTTACAATAGCTTCTCCTGTTATTG...
    >BHYC-4691
    GAAGTGAACTGTGTTGCTTGTGGGTTTCTTGCTGCTCTTGCTGTCACTGCTTCTCCCGTAATCGCTG...
    etc...

### Read files

You will need to provide the `hybpiper.nf` pipeline with either:

a) A directory of paired-end forwards and reverse reads (and, optionally, a file of single-end reads) for each sample; or

b) A directory of single-end reads.

For the read files to be recognised by the pipeline, they should be named according to the default convention:

    *_R1.fastq 
    *_R2.fastq
    *_single.fastq (optional; will be used if running with the flag `--paired_and_single`)

OR

    *_R1.fq 
    *_R2.fq
    *_single.fq (optional; will be used if running with the flag `--paired_and_single`)

**NOTE:**

- It’s fine if there’s text after the `R1`/`R1` and before the `.fastq`/`.fq`, as long as it's the same for both read files.
- It’s fine if the input files are gzipped (i.e. suffix `.gz`). 
- You can provide a folder of single-end reads only (use the flag `--single_only` if you do). Your read files should be named `*_single.fastq` (or `.fastq.gz`, `.fq`, `.fq.gz`).  
- You can specify a custom pattern used for read file matching via the parameters `--read_pairs_pattern <pattern>` or `--single_pattern <pattern>`.
- If your samples have been run across multiple lanes, you'll likely want to combine read files for each sample before processing. The pipeline can do this for you - see [here][22] for details.

## Running on Linux

Please see the Wiki entry [Running on Linux][3].

## Running on a Mac (macOS)

Please see the Wiki entry [Running on a Mac][1].

**NOTE:** Macs using the new Apple M1 chip are not yet supported.

## Running on a PC (Windows)

Please see the Wiki entry [Running on a PC][20].

## Nextflow pipeline options and parameters


Example run command:

    nextflow run hybpiper.nf -c hybpiper.config -entry assemble -profile standard_singularity --illumina_reads_directory reads_for_hybpiper --targetfile_dna Angiosperms353_targetSequences.fasta

```   
\Mandatory arguments:

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
                              tab-delimited file (one <gene>\t<taxon_name> per line)
                              or a single taxon name

  --exclude                   Do not use any sequence with the specified taxon name
                              string in Exonerate searches. Sequenced from this
                              taxon will still be used for read sorting

  --timeout_extract_contigs   Kill long-running processes if they take longer than
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
```

Please see the Wiki entry [Additional pipeline features and details][5] for further explanation of the parameters above, and general pipeline functionality.

## Output folders and files

If you're just after the unaligned `.fasta` files for each of your target genes (not including putative paralogs), the two main output folders of interest are probably:

- `07_sequences_dna`
- `08_sequences_aa`

If you need the per-sample output folders produced by a standard HybPiper run, these can be found in folder:

- `04_processed_gene_directories`

For a full explanation of output folders and files, please see the Wiki entry [Output folders and files][6].

## General information

For details on adapting the pipeline to run on local and HPC computing resources, see [here][18].

## Issues still to deal with (WIP)

Please see the Wiki entry [Issues][4].

## Changelog

*11 November 2025*

- Bugfix: syntax 'optional true' changed to ', optional: true' in `hybpiper.nf` script.
- Updated the `hybpiper-nf` script to version 1.1.1

*30 April 2025*

- Updated the HybPiper version in the Singularity container `hybpiper-paragone` to version 2.3.2, and updated the `hybpiper-nf` script to version 1.1.0.
- Updated the command line options for the entry point 'assemble' to reflect HybPiper version 2.3.2.

*31 May 2024*

- Updated the HybPiper version in the Singularity container `hybpiper-paragone` to version 2.1.7, and updated the `hybpiper-nf` script to version 1.0.4
- Fixed an issue when calling `hybpiper retrieve_sequences` with an amino-acid target file.
- Updated the flag `--use_diamond` to `--diamond` to match stand-alone HybPiper.
- Updated the flag `--run_intronerate` to `--no_intronerate` to reflect changes made in HybPiper >= 2.1.6


*03 July 2023*

- Updated the HybPiper version in the Singularity container `hybpiper-paragone` to version 2.1.5, and updated the `hybpiper-nf` script to version 1.0.3
- Added the new options `exonerate_hit_sliding_window_size` and `exonerate_hit_sliding_window_thresh` to the `hybpiper assemble` options.

*15 May 2023*

- Bugfix: BBmap.sh memory for Hybpiper chimera test was set to a default of 1 Mb. Changed to 1000 Mb.
- Update version in hybpiper-nf script to version 1.0.2

*14 April 2023*

- Updated the HybPiper version in the Singularity container `hybpiper-paragone` to version 2.1.3, and updated the `hybpiper-nf` script to version 1.0.1
- Bugfix: get target file basename for HybPiper commands (fixes error when using a target file that isn't in the current working directory)

*01 February 2023*

- Added a `conda` and `conda_slurm` profile. This allows the pipeline to be run using conda packages rather than the Singularity container. The corresponding conda environment is created in the nextflow `work` directory.
- Updated the `hybpiper.nf` script to support all native HybPiper 2 parameters.
- Added a script version number; view by using the flag `--version`.
- Split the main HybPiper `assemble` pipeline, and the HybPiper commands `check_targetfile` and the `fix_targetfile`; these are now separate entry points to the `hybpiper.nf` script, accessible by using the parameter `-entry`, e.g. `-entry assemble`.
- Added separate help docs for each entry point, e.g. `-entry assemble --help`.

*28 November 2022*

- Change repository name from `HybPiper-RBGV` to `hybpiper-nf`.
- Update the required container name from `hybpiper-yang-and-smith-rbgv.sif` to `hybpiper-paragone.sif`.
- Update containerised HybPiper version and corresponding Nextflow scripts to HybPiper version 2

*09 November 2021* 

- Nextflow script updated to DSL2. **NOTE:** Nextflow version >= 21.04.1 is now required!
- Singularity *.def file updated to use Ubuntu 20.04, and to install Muscle and FastTreeMP.

[1]:https://github.com/chrisjackson-pellicle/hybpiper-nf/wiki/Running-on-a-Mac-(macOS)-with-Vagrant "Link to Running on a Mac Wiki"
[3]:https://github.com/chrisjackson-pellicle/hybpiper-nf/wiki/Running-on-Linux "Link to Running on Linux Wiki"
[4]:https://github.com/chrisjackson-pellicle/hybpiper-nf/wiki/Issues-(WIP) "Link to Issues Wiki"
[5]:https://github.com/chrisjackson-pellicle/hybpiper-nf/wiki/Additional-pipeline-features-and-details "Link to Additional pipeline features and details Wiki"
[6]:https://github.com/chrisjackson-pellicle/hybpiper-nf/wiki/Output-folders-and-files "Link to Output folders and files Wiki"
[7]:http://www.usadellab.org/cms/?page=trimmomatic "Link to Trimmomatic website"
[8]:https://www.ncbi.nlm.nih.gov/books/NBK279690/ "Link to BLAST command line documentation"
[9]:https://biopython.org/ "Link to BioPython website"
[10]:https://github.com/ablab/spades "Link to SPAdes assembler website"
[11]:http://bio-bwa.sourceforge.net/ "Link to BWA mapper website"
[12]:https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/ "Link to bbmap documentation"
[13]:https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate "Link to Exonerate website"
[14]:http://www.htslib.org/ "Link to htslib website for Samtools"
[15]:https://sylabs.io/docs/ "Link to Singularity website"
[16]:https://www.nextflow.io/ "Link to Nextflow website"
[17]:https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4948903/ "Link to HybPiper manuscript"
[18]:https://github.com/chrisjackson-pellicle/hybpiper-nf/wiki/Additional-pipeline-features-and-details#managing-computing-resources "Link to managing computing resources"
[19]:https://github.com/mossmatters/HybPiper/wiki#12-target-file "Link to HybPiper Wiki target file details"
[20]:https://github.com/chrisjackson-pellicle/hybpiper-nf/wiki/Running-on-a-PC-(Windows)-with-Vagrant "Link to Running on a PC Wiki"
[21]:https://slurm.schedmd.com/overview.html "Link to SLURM website"
[22]: https://github.com/chrisjackson-pellicle/hybpiper-nf/wiki/Additional-pipeline-features-and-details#combining-read-files-for-samples-run-across-multiple-lanes "Link to Wiki section combining-read-files-for-samples-run-across-multiple-lane"
