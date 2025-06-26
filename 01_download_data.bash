#!/bin/bash
# Laura Dean
# 17/6/25
# script written for running on the UoN HPC Ada

#SBATCH --job-name=DepMap_download
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10g
#SBATCH --time=10:00:00
#SBATCH --output=/gpfs01/home/mbzlld/code_and_scripts/slurm_out_scripts/slurm-%x-%j.out

# set variables
wkdir=/gpfs01/home/mbzlld/data/rembrandt
cell_lines=$wkdir/cell_lines.txt # tab delimited file, with header line, containing cell line names followed by run accessions
cell_lines=$wkdir/SF188_accessions.txt
# move to working directory
cd $wkdir


# Download the rembrandt data with wget
#wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE108nnn/GSE108476/matrix/GSE108476-GPL570_series_matrix.txt.gz
# extract the data
#gunzip GSE108476-GPL570_series_matrix.txt.gz
# remove this dataset as it isn't useful for our purposes
#rm GSE108476-GPL570_series_matrix.txt
# This downloaded the expression data fine but I then realised that the rembrandt dataset
# is array based and information is only present at the gene level, so we cannot use this
# database to search for splice variants. I will try with the DepMap database as this is
# RNAseq and also has glioblastoma samples.


# Download the DepMap raw data
# an example for a single cell line
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/009/SRR8615649/SRR8615649_1.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR861/009/SRR8615649/SRR8615649_2.fastq.gz

# for multiple cell lines from DepMap
# create an array of accession numbers from input file that contains cell line names followed by run accessions
accessions=($(tail -n +2 $cell_lines | cut -f2))

# change between wget commands below depending on whether you data is paired or single end
for srr in ${accessions[@]}; do
  pre=${srr:0:6}
  end="00${srr: -1}"
  #wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$pre/$end/$srr/${srr}_1.fastq.gz
  #wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$pre/$end/$srr/${srr}_2.fastq.gz
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$pre/$end/$srr/${srr}.fastq.gz
done

