import matplotlib
matplotlib.use('Agg')
import os,datetime
import os.path
import click
import MIST
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 
import sys   
sys.setrecursionlimit(100000)
import faulthandler;faulthandler.enable()
@click.group()
def cli():
    """MIST - Metagenomic Intra-Species Typing
    """
    pass



@cli.command()
@click.option('-p','--threads',type=int,help="Number of threads for Bowtie2 (default: 8)")
@click.option('-1','--pair_1',type=click.Path(exists=True),help="input fq file with #1 mate, paired with pair_1 file")
@click.option('-2','--pair_2',type=click.Path(exists=True),help="input fq file with #2 mate, paired with pair_2 file")
@click.option('-d','--database',type=click.Path(exists=True),help="input bowtie2 index file for the pan-genome sequences")
@click.option('-o','--output',type=click.Path(exists=False, writable=True),help="output folder which contains: 1) read counts for each pathogen species (_MIST_species_count.txt); 2) reads specific to each pathogen species (_MIST.*.fq).")
def species(pair_1,pair_2,database,output,threads=8):
    """- conduct species-level typing"""
    MIST.species(threads,pair_1,pair_2,database,output)

@cli.command()
@click.option('-i','--refdir',type=click.Path(exists=True),help=" input folder of reference genomes")
@click.option('-o','--output',type=click.Path(exists=False, writable=True),help="output folders containing index files")
def index(refdir,output):
    """- build index files for reference genomes"""
    MIST.index(refdir,output)

@cli.command()
@click.option('-i','--refdir',type=click.Path(exists=True),help="input folder of reference genomes")
@click.option('-o','--output',type=click.Path(exists=False, writable=True),help="matrix file of the clustered reference genomes")
@click.option('-t','--threads',type=int,help="Number of threads for ANI")
@click.option('-s','--cutoff',type=str,help=" list of similarity thresholds (between 0 and 1); separated by comma (e.g. 0.98,0.99)")
def cluster(refdir,output,threads=1,cutoff=[0.98,0.99,0.999]):
    """- cluster reference genomes by ANI"""
    MIST.cluster(refdir,output,threads,cutoff)

@cli.command()
@click.option('-p','--threads',type=int,help="Number of threads for Bowtie2")
@click.option('-i','--indexpath',type=click.Path(exists=True),help="input folder of reference index files")
@click.option('-U','--single_end',type=click.Path(exists=True),help="input single-end fq file")
@click.option('-1','--pair_1',type=click.Path(exists=True),help="input fq file with #1 mate, paired with pair_1 file")
@click.option('-2','--pair_2',type=click.Path(exists=True),help="input fq file with #2 mate, paired with pair_2 file")
@click.option('-l','--read_length',type=click.FLOAT,help="read length")
@click.option('-o','--output',type=click.Path(exists=False, writable=True),help="output folder for mismatch matrix file and alignment output files")
def map(indexpath,read_length,output,threads=8,single_end=None,pair_1=None,pair_2=None):
    """- map mNGS reads against reference genomes"""
    MIST.map(indexpath,read_length,output,threads,single_end,pair_1,pair_2)

@cli.command()
@click.option('-c','--cluster_output',type=click.Path(exists=True),help="input file; the matrix of the clustered reference genomes")
@click.option('-m','--mismatch_matrix_output',type=click.Path(exists=True),help="input file; the mismatch matrix of the alignments")
@click.option('-l','--read_length',type=click.FLOAT,help="read length")
#@click.option('-a','--a_value',type=click.FLOAT,help="a value")
@click.option('-o','--output',type=click.Path(exists=False, writable=True),help="output folder for result files")
def subspecies(cluster_output,mismatch_matrix_output,read_length,output):
    """- estimate strain-level composition"""
#    start = time.time()
    MIST.subspecies(cluster_output,mismatch_matrix_output,read_length,output)
#    end = time.time()
#    print("subspecies module run time:",end-start)
if __name__=="__main__":
    cli()

