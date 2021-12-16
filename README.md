
# MIST (Metagenomic Intra-Species Typing)

## Introduction
MIST is a metagenomic intra-species typing technique that was developed primarily for clinical specimens with low pathogen loads. Hopefully, it will aid in strain-level diagnoses of bacterial infections as well as public health epidemiology and surveillance. Its algorithm contains the following three characteristics. 
* Based on average nucleotide identity (ANI), reference genomes are clustered into hierarchical levels to resolve the ambiguous definition of “strain”.
* Maximum likelihood estimation is conducted upon the reads’ mismatch values to infer the compositional abundance.
* Read ambiguity is used to infer the abundance uncertainty, and the similarity to reference genomes is used to predict the presence of novel strains.
MIST contains four modules: Species, Cluster, Index, Strain. Its workflow is depicted in the figure below. MIST depends on the counts of reads that are uniquely matched to the pan-genomes of each pathogen species for __species-level typing__ (panel A). Next, MIST prepares a hierarchical database of reference genomes based on ANI grouping for __strain-level typing__ (panel B). By matching reads to all of the species' reference genomes, the scores for each alignment are transformed to posterior probabilities that indicate the likelihood of the sequence read of being allocated to each cluster. The probability matrix is then employed using the maximum likelihood estimation (MLE) to infer the abundance of each cluster.

<p align="center"><img src="https://github.com/pandafengye/MIST/blob/master/Pipeline.png" alt="workflow_small"  width="800">

__Note: Red box, the module "Species"; yellow box, the modules "Index" and "Cluster"; purple box, the module "Strain".__

## Requirements
### Software
* Linux system
* [Python = 3.6](https://www.python.org)
* [GCC >= 4.8](https://gcc.gnu.org)
* [Bowtie 2](https://github.com/BenLangmead/bowtie2)
* [FastANI](https://github.com/ParBLiSS/FastANI)
* Python modules: `networkx`, `pandas`, `matplotlib`, `numpy`, `scipy`, `scikit-learn`, `joblib`, `click`.

### Pre-built database
* Pre-built-pangenome: This folder is used in the “species” module, which include bowtie-indexed pan-genomes of 14 bacterial species (listed below).
* Pre-built bowtie-indexed reference genomes: For each species, their complete genomes deposited from NCBI Genbank database are downloaded and are bowtie-indexed with the module “index”. These pre-built bowtie index files are used in combination with the pre-built clustering files in the module “strain”.
   * [Acinetobacter_baumannii](http://bacdb.cn/Acinetobacter_baumannii_MIST_index.tgz)
   * [Campylobacter_jejuni](http://bacdb.cn/Campylobacter_jejuni_MIST_index.tgz)
   * [Clostriodioides_difficile](http://bacdb.cn/Clostriodioides_difficile_MIST_index.tgz)
   * [Enterococcus_faecalis](http://bacdb.cn/Enterococcus_faecalis_MIST_index.tgz)
   * [Enterococcus_faecium](http://bacdb.cn/Enterococcus_faecium_MIST_index.tgz)
   * [Escherichia_coli](http://bacdb.cn/Escherichia_coli_MIST_index.tgz)
   * [Haemophilus_influenzae](http://bacdb.cn/Haemophilus_influenzae_MIST_index.tgz)
   * [Klebsiella_pneumoniae](http://bacdb.cn/Klebsiella_pneumoniae_MIST_index.tgz)
   * [Legionella_pneumophila](http://bacdb.cn/Legionella_pneumophila_MIST_index.tgz)
   * [Listeria_monocytogenes](http://bacdb.cn/Listeria_monocytogenes_MIST_index.tgz)
   * [Mycobacterium_tuberculosis](http://bacdb.cn/Mycobacterium_tuberculosis_MIST_index.tgz)
   * [Salmonella_enterica](http://bacdb.cn/Salmonella_enterica_MIST_index.tgz)
   * [Staphylococcus_aureus](http://bacdb.cn/Staphylococcus_aureus_MIST_index.tgz)
   * [Streptococcus_pneumoniae](http://bacdb.cn/Streptococcus_pneumoniae_MIST_index.tgz)

* Pre-built-clustering files: This folder contains the clustering files (obtained from the module “cluster”) for the 14 pathogens and is used in the module “strain”. Same as above, these clustering files can only be used in combination with the bowtie-indexed reference genomes.

__Note:__ In additional to the pre-built database above, you can customize your own database by the following steps:
   * Downloading reference genomes of a certain species in FASTA format.
   * Cluster these reference genomes by the “cluster” module.
   * Bowtie-index these reference genomes by the “index” module.


## Installation
  ```bash
  $ conda create -n MIST -c conda-forge -c bioconda python=3.6 fastANI bowtie2
  $ conda activate MIST
  $ git clone https://github.com/pandafengye/MIST.git
  $ cd MIST
  $ pip install -r requirements.txt --default-timeout=1000 # Install related python dependencies
  $ python MIST.py # Test install
  ```
## Usage
MIST contains four modules: __index__, __cluster__, __species__ and __strain__.

## MIST-index
  This module functions to index the reference genomes with Bowtie2 indexer (bowtie2-build). Once the reference genomes are indexed, users will not need to re-index the
genome before each analysis of metagenomics datasets.
### Command
  ```bash
  $ python MIST.py index --refdir Example_Dir/input/ref_dir/ --output Example_Dir/output/
  ```
### Options: 
    -i, --refdir PATH 
         Path to the reference genome folder; All reference genomes should be in FASTA format and put in the same folder; each file represents one reference genome, with a *.fa prefix.
    -o, --output PATH
       Output folder saving the index files for reference genomes. The base name of the index is the same as the reference genome.
 ## MIST-cluster
  This module functions to assign reference genomes into clusters at user-defined levels. MIST calls FastANI program to calculate ANI for estimation of pairwise genetic distance, based on which the reference genomes are divided into clusters. Same as the Index module, once the clusters are established, users are not required to run this module before each independent job.

### Command
  ```bash
  $ python MIST.py cluster --threads 8 --refdir Example_Dir/input/ref_dir/ --cutoff 0.98,0.99,0.999 --output Example_Dir/output/
  ```
### Options:
    -t, --threads INTEGER  Number of threads for ANI
    -i, --refdir PATH     
    input folder of reference genomes
    -o, --output PATH    
    matrix file of the clustered reference genomes
    -s, --cutoff TEXT    
    list of similarity thresholds (between 0 and 1); separated by comma (e.g. 0.98,0.99)

## MIST-species
This module functions to perform species-level typing. MIST calls Bowtie2 to map the user’s mNGS reads (in .fastq format) against the pan-genomes of each bacterial species
and estimate the abundance by counting the reads mapped to each species. The species-specific reads are extracted from the resulting SAM file for the downstream strain-level
typing.
### Command
```bash
$  python MIST.py species --threads 8 --pair_1 Example_Dir/input/read/example_data1.1.fq --pair_2 Example_Dir/input/read/example_data1.2.fq --database Pre-built-pangenome/ --output Example_Dir/output/
```
### Options:
  	-p, --threads INTEGER      
	Number of threads for Bowtie2 (default: 8)
  	-1, --pair_1 PATH         
	input fq file with #1 mate, paired with pair_2 file
  	-2, --pair_2 PATH
  	input fq file with #2 mate, paired with pair_1 file
  	-d, --database PATH
	input bowtie2 index file for the pan-genome sequences
  	-o, --output PATH        
	output folder which contains: 1) read counts for each pathogen species (_MIST_species_count.txt); 2) reads specific to each pathogen species (_MIST.*.fq).
### Tips:
  The pre-built pan-genome index file is available at http://bacdb.cn/Pre-built-pangenome.tgz. 
  For the reads specific to each pathogen species (\_MIST.\*.fq), 0.1x sequencing coverage of bacterial genome (e.g. 5000 100-bp reads for a 5-Mb bacterial genome) is usually sufficient for MIST to do strain-level typing. Too many reads (e.g., > 50000 reads) for the subsequent mapping and maximum likelihood estimation would otherwise cause long running time. Users can extract a subset (5000) of reads with the command such as "head –n 20000 _MIST.*.fq > input.fq".


## MIST-strain
  This module functions to map metagenomic sequences against reference genomes using Bowtie2, and to measure the relative abundance of each cluster in the metagenomics dataset, along with similarity and reliability assessment.

### Command
  ```bash
  $ python MIST.py strain --threads 8 --indexpath Example_Dir/output/_MIST_index/ --cluster_output Example_Dir/output/_MIST_ref_cluster.csv --pair_1 Example_Dir/input/read/test.1.fq --pair_2 Example_Dir/input/read/test.2.fq --read_length 200 --output Example_Dir/output/
  ```
### Options:
      -p, --threads INTEGER      
    Number of threads for Bowtie2 (default: 8)
      -c, --cluster_output PATH  
     input file; the matrix of the clustered reference genomes
      -i, --indexpath PATH       
    input folder of index files for reference genomes; produced by MIST-index module.
      -U, --single_end PATH     
    input single-end fq file;can be the output produced by MIST-species module.
      -1, --pair_1 PATH         
    input fq file with #1 mate, paired with pair_2 file
      -2, --pair_2 PATH         
    input fq file with #2 mate, paired with pair_1 file; either choose the paired input or the single-end input.
      -l, --read_length FLOAT   
    read length
      -g, --genome_size INTEGER  
    genome size (optional)
      -o, --output PATH        
    output folder for mismatch matrix file and alignment output files. A folder _MIST_map_alignment, which contains the mapped .sam files corresponding to each reference genome; a file _MIST_map_Mismatch_matrix.csv, which contains the number of mismatches derived from each read mapping against each reference genome.
