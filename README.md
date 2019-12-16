# MIST （Metagenomic Intra-Species Typing）
  MIST is a mmetagenomic intra-species typing tool which is specifically designed for clinical specimens of low pathogen load. Its algorithm has the following three features.1) Based on average nucleotide identity (ANI), reference genomes are clustered into hierarchical levels to resolve the ambiguous definition of “strain”; 2) Maximum likelihood estimation is conducted upon the reads’ mismatch values to infer the compositional abundance. 3) Read ambiguity is used to quantify the abundance uncertainty along with the similarity of reference genomes to specimen’s constituents. Hopefully, it benefits strain-level diagnostics as well as public health epidemiology and surveillance.
# Prerequisites
  MIST pipeline was built and tested in the Intel architectures (x86_64) running Linux/Unix      environment. It may work under other operating systems and have not been tested.
Bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and FastANI (https://github.com/ParBLiSS/FastANI) are required to be preinstalled. The locations of their binaries should be included in the environment path ($PATH).
# Download and install 
  ```bash
  $ git clone https://github.com/pandafengye/MIST.git
  $ cd MIST
  $ python MIST_cmd.py
  ```
# Usage
MIST contains four modules: index, cluster, map, measure.
## MIST-Index
  This module functions to index the reference genomes with Bowtie2 indexer (bowtie2-build). Once the reference genomes are indexed, users do not need to re-index the genome before each analysis of metagenomics datasets.
### Command
  ```bash
  $ python MIST_cmd.py index -i test/sampledata/ref_dir/ -o test/result/
  ```
### Options: 
    -i, --refdir PATH 
         Path to the reference genome folder; All reference genomes should be in FASTA format and put in the same folder; each file     represents one reference genome, with a *.fa prefix. See the example folder ‘test/sampledata/ref_dir’.
    -o, --output PATH
       Output folder saving the index files for reference genomes. The base name of the index is the same as the reference genome.
 ## MIST-cluster
  This module functions to assign reference genomes into clusters at user-defined levels. MIST calls FastANI program to calculate ANI for estimation of pairwise genetic distance, based on which the reference genomes are divided into clusters. Same as the Index module, once the clusters are established, users are not required to run this module before each independent job.

### Command
  ```bash
  $ python MIST_cmd.py cluster -i test/sampledata/ref_dir/ -s 0.95,0.98,0.99 -o test/result/
  ```
### Options:
    -i, --refdir PATH     
    input folder of reference genomes
    -o, --output PATH    
    matrix file of the clustered reference genomes
    -s, --cutoff TEXT    
    list of similarity thresholds (between 0 and 1); separated by comma (e.g. 0.98,0.99)

## MIST-map
  This module functions to map metagenomic sequences against reference genomes using Bowtie2. 
### Command
  ```bash
  $ python MIST_cmd.py map -p 8 -i test/result/index/ -1 test/sampledata/read/test.1.fq -2 test/sampledata/read/test.2.fq -l 100 -o test/result/ 
  ```
### Options:
      -p, --thread INTEGER      
    Number of threads for Bowtie2 (default: 8)
      -i, --indexpath PATH       
    input folder of index files for reference genomes; produced by MIST-index module.
      -U, --single_end PATH     
    input single-end fq file
      -1, --pair_1 PATH         
    input fq file with #1 mate, paired with pair_2 file
      -2, --pair_2 PATH         
    input fq file with #2 mate, paired with pair_1 file; either choose the paired input or the single-end input.
      -l, --read_length FLOAT   
    read length
      -o, --output PATH        
    output folder for mismatch matrix file and alignment output files. A folder _MIST_map_alignment, which contains the mapped .sam files corresponding to each reference genome; a file _MIST_map_Mismatch_matrix.csv, which contains the number of mismatches derived from each read mapping against each reference genome. See example fold ‘test/result/’.
## MIST- measure
  This module functions to measure the relative abundance of each cluster in the metagenomics dataset, along with similarity and reliability assessment.
### Command
```bash
$ python MIST_cmd.py measure -c test/result/_MIST_ref_cluster.csv -m test/result/_MIST_map_Mismatch_matrix.csv -l 100 -o test/result/
```
### Options:
    -c, --cluster_output PATH  
    input file; the matrix of the clustered reference genomes
    -m, --mismatch_matrix_output PATH 
    input file; the mismatch matrix of the alignments
    -l, --read_length FLOAT   
    read length
    -o, --output PATH        
    output folder in which the result files correspond to each level of cluster. Information of the estimated abundance, 95% CI, P value and similarity are given for each cluster in each file. See example fold ‘test/result/’.




