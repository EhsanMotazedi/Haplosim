HAPLOSIM, version 1.8, Ehsan Motazedi (ehsan.motazedi@gmail.com), last modified: 09 12 2017

Introduction
============

**Haplosim**\(*haplosim.py*\) is a simulation pipeline to evaluate the performance of single individual haplotyping algorithms that make use of (polyploid) next generation sequencing (NGS) reads. Having specified a reference DNA sequence in fasta format, mutations are induced at single nucleotide level, according to the chosen stochastic model, by software HaploGenerator to simulate individual genomes with the desired ploidy level. According to the chosen technology, NGS reads are generated for Illumina (GA1, GA2, MiSeq, HiSeq 2000, HiSeq 2500) or PacBio (CCS and CLR) in silico, and haplotypes are reconstructed from those reads by three algorithms: HapCompass (Aguiar and Istrail 2013), HapTree (Berger *et al.* 2014) and SDhaP (Das and Bansal 2015). The accuracy of the estimated haplotypes will be assessed against the original haplotypes using several measures, e.g. Vector Error Rate (VER), Phasing Accuracy Rate (PAR), Proportion of missing variants, etc. As an option, one could also obtain compuation timeof the CPU and the memory consumption of each algorithm for each simulated genome. The quality evaluations for each algorithm are reported in the algorithm's DAT file, and the time/memory analysis is reported in separate "timem" DAT file. 

The pipeline is compatible with Python 2.7 or later versions, but NOT with Python 3, and makes use of several BASH commands. It is therefore designed for POSIX-oriented platforms, e.g. UNIX and Mac-OS.

Requirements:
=============

To run *Haplosim*, you need the following input file(s):

1. A reference genome in fasta format from which artificially modified genomes are produced. The name of the target contig, length and coordinate of the modified genomes could be specified as input parameters. In the random mode, genomic regions with the desired length are selected at random sites to generate the mutants. In this case, the number of randomly selected regions and the number of mutants for each may be specified. 

Besides, you need the following softwares (and their requirements) installed on your system:

	1. ART to simulate Illumina reads (NOTE: *Haplosim* currently supports only **ChocolateCherryCake** version!): 

	http://www.niehs.nih.gov/research/resources/software/biostatistics/art/

	At least one of the following aligners:

  		2. Blasr (for PacBio)

  		https://github.com/PacificBiosciences/blasr

  		3. Bowtie 2 (for Illumina):

  		http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

  		4. BWA (Burrows-Wheeler Aligner) (for both Illumina and PacBio):

  		http://bio-bwa.sourceforge.net/

	5. ExtractHAIRS from HapCUT to generate the fragment matrix from VCF file:

	https://github.com/vibansal/hapcut

	6. FreeBayes to call the variants from the alignment:

	https://github.com/ekg/freebayes

	7. HapCompass_v0.8.2 to run HapCompass algorithm:

	http://www.brown.edu/Research/Istrail_Lab/hapcompass.php

	8. HapTree_v0.1 to run HapTree algorithm (wrapped by the python script HapTree_multiblock.py, see HapTree_multiblock --help):

	http://groups.csail.mit.edu/cb/haptree/

	9. PBSIM to simulate PacBio reads:

	https://code.google.com/p/pbsim/

	10. Picardtools (AddOrReplaceReadGroups):

	http://broadinstitute.github.io/picard/

	11. Samtools for indexing, sorting and duplicate removal (assumes version 1.4):

	https://github.com/samtools/samtools

	12. SDhaP to run SDhaP algorithm (wrapped by the bash script hap2, see hap2 --help):

	http://sourceforge.net/projects/sdhap/

Java must be also installed on the system, as it is necessary for running HapCompass and Picardtools. 

Format of the output files:
================================

1. Three DAT files with the following fields, each line reports the values for one simulation, i.e. region- mutant:

	**CORR_ALL_DOSAGES**: Allelic correlation score between between the estimated and original haplotypes, comparing all the variants common between the original and estimates. This measure is equal to the probability that alleles at the same position are the same between the two haplotypes.

	**PSR_ALL_DOSAGES**: Pairwise-phasing Accuracy Rate of the estimated haplotypes, comparing all the variants common between the original and estimates. This is the proportion of correct pairwise-phasings among the phasings of all possible pairs of SNPs. 

	**VER**: Vector Error Rate (Berger et al. 2014) of the estimated haplotypes, including the common variants whose dosage is correctly estimated, normalized by the number of these variants and the ploidy.

	**PSR_CORRECT_DOSAGES**: like PSR_ALL_DOSAGES, excluding variants with wrongly estimated dosages.

	**CORR_CORRECT_DOSAGES**: like CORR_ALL_DOSAGES, excluding variants with wrongly estimated dosages.

	**MISSING_ESTIMATE**: Proportion of original SNPs missing in the estimate.

	**SPURIOUS_ESTIMATE**: Ratio of spurious (false) variants in the estimate to the original variants. 

	**WRONG_DOSAGES_COMMON_SNPS**: Proportion of the common SNPs with wrong dosage in the estimate.

	**PPV_ESTIMATE**: Positive predictive value of the output SNPs.

	**NUM_GAPS**: Number of gaps induced during the estimation of the haplotyoes, i.e. number of output haplotype blocks minus one, normalized by the number of SNPs.

	**CORRECT_HOMOLOGUE_RATE**: The proportion of the correctly estimated homologues in the simulated population.

	**ID_REGION**: ID of the reference genomic region used to simulate the mutant genome.

	**ID_MUTATION**: Individual IDs of the simulated genomes. This was NOT given in version 1.0.


Also, in the --report-aln mode:

2. fasta reference, fastq reads, aligned BAM ans SAM and VCF files for each simulation, names as my{ref, reads, bam, sam, vcf}xy.{fa, fastq, bam, sam, vcf} for selected reference region number x and mutant numebr y.: This file should contain the variants and genotypes for a single individual that is to be phased. Please make sure that the VCF file conforms to the standard VCF format, since HapCUT does not check this. 

3. Log file containing all the generated and estimated haplotypes. Haplotypes only include the variant sites according to definition. The variant position (the same as in the VCF file), ref and alternative alleles are specified in the haplotypes. Each line represents a variant and its allelic dosage in the haplotypes. 

In the verbose mode, the number of failures for each algorithm, as well as the error and warning messages and the total CPU time are written to the log file. If --report-aln is off, only problematic haplotypes are reported in the log.

Also, in the verbose mode:

4. A separate log file for the read generator, including the output, error and warning messages from ART and PBSIM. 

Running the program:
=====================

Before Running the program, the path to the required programs must be set in the OS environment, e.g by exporting following variables:

	ARTSPath, blasrPath, bowtiePath, bwaPath, FreebayesPath, HapCompassPath, HapCutPath, HapTreePath, PBSIMPath, PicardtoolsPath, samtoolsPath, SDhaPPath. 

The default is to $HOME directory and the default sub directory containing the executable files. Note that the default folder name might be different from the actual folder on your system. It is also possible to add the path to the executable files to $PATH in the OS environment, hence not needing the variables above (except for *Picardtools* and *HapCompass* which are run as java scripts and whose correct paths must be exported to *Haplosim*).
    
For parallel execution of *HapTree\_multiblock*, NCORES variable must be set in the OS environment to the desired number of cores.    

Input parameters:
=====================

For the input parameters, see the help of the software:

./haplosim.py --help, -h

***Examples:***


Simulate tetraploid genomes (-p 4) from chromosome 5 of PGSC\_DM\_v4.03 draft potato genome (Sharma *et al.* 2011), using a lognormal model with mean and variance of the log distance between neighboring SNPs set to 3 and 1, respectively. 25 regions of 20kb length (-l 20e03) are randomly chosen from this reference from each 20 genomes are generated by random mutation (-r 25 20) according to the lognormal model. Paired-End Reads are generated by MiSeq technology with insert size 600, with a 15x coverage per homologue. The proportions of simplex, duplex, triplex and tetraplex alleles are specified by --dosage. Results will be stored in potato_15_MS_600_lognormal_20e03_HapCompass.dat, potato_15_MS_600_lognormal_20e03_HapTree.dat, potato_15_MS_600_lognormal_20e03_SDhaP.dat and 
potato_15_MS_600_lognormal_20e03.log files.

    haplosim.py -x 15 $HOME/PGSC_v4.03_resources/PGSC_DM_v4.03_CHR05.fasta potato_15_MS_600_lognormal_20e03 MS PE --insert-size 600 -p 4 -l 20e03 -r 25 20 --lognormal 3 1 --dosage 0.5 0.23 0.14 0.13

Change the coverage to 10x per homologue, stochastic model to mixture poisson with main mutation rate 0.04, technology to PacBio long reads. Aligner is set to bwa-mem (default aligner is blasr for PacBio) and the time out to 1200 seconds (default 500s). --no-missing option makes the software keep all of the estimated haplotypes reported. Dosage distribution will be set to default, i.e. equal probabilities for all of the alternative allele counts simplex, douplex, triplex and tetraplex.

    haplosim.py -x 10 $HOME/PGSC_v4.03_resources/PGSC_DM_v4.03_CHR05.fasta potato_10_PacBio_poisson_20e03 CLR -p 4 -l 20e03 -r 15 20 -m 0.04 -t 1200 --no-missing --bwa-mem

A limit could be also set on the total number of generated haplotypes from each reference, using --maxhap option. In that case, random haplotypes are generated equal to the number passed to --maxhap, and the homologues of each genome are selected by sampling with replacement from this total set.


**HaploGenerator** (*haplogenerator.py*) could be used independently to generate random genomes from a reference, by inducing mutations and indels, with several stochastic models. To produce allopolyploid genomes, sub-genomes could be generated from the same reference genome independently by HaploGenerator, and then combined by **allowrapper**. The output will be the fasta file (and its index) for each homologue, as well as a text file containing the generated haplotypes. 

***Examples:***

Generate hexaploid genome from tetraploid and diploid sub-genomes, with the names genHexa\_hap1, ..., genHexa\_hap6, genHexa\_varianthaplos.txt:

    haplogenerator.py -f $HOME/myref.fa -o genA --model poisson -s "[0.04,0,0]" -m "{'A':'C','G':'A','C':'T','T':'G'}" -p 4 

    haplogenerator.py -f $HOME/myref.fa -o genB --model lognormal -s "[3,9,0]" -m "{'A':'T','G':'C','C':'A','T':'GAC'}" -p 2 --sdlog "[1,1,0]"

    allowrapper -p 6 -f $HOME/myref.fa genA genB genHexa 

    
Wrapper script is provided for HapTree \(*HapTree_multiblock.py*\) to allow multiple haplotype block. Wrapper is provided for SDhaP \(*hap2*\) to reformat its output to be comparable to the outputs of other haplotyping software's and HaploGenerator. 

**HapCompare** \(*hapcompare.py*\) could be used independently to compare haplotype files in the format of the outputs of *HaploGenerator*, *HapTree\_multiblock*, *hap2* and *HapCompass*.


Download:
=====================
The software is available at gitlab, Wageningen UR:

www.bif.wur.nl
