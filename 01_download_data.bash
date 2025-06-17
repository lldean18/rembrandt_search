#!/bin/bash
# Laura Dean
# 17/6/25
# script written for running on the UoN HPC Ada

#SBATCH --job-name=rembrandt_download
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10g
#SBATCH --time=10:00:00
#SBATCH --output=/gpfs01/home/mbzlld/code_and_scripts/slurm_out_scripts/slurm-%x-%j.out

# set variables
wkdir=/gpfs01/home/mbzlld/data/rembrandt

# move to working directory
cd $wkdir

# Download the rembrandt data with wget
wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE108474&format=file


