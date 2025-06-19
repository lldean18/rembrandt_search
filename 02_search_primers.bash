#!/bin/bash
# Laura Dean
# 19/6/25
# script written for running on the UoN HPC Ada

#SBATCH --job-name=bbmap
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10g
#SBATCH --time=4:00:00
#SBATCH --output=/gpfs01/home/mbzlld/code_and_scripts/slurm_out_scripts/slurm-%x-%j.out

# set variables
wkdir=/gpfs01/home/mbzlld/data/rembrandt
reads1=$wkdir/SRR8615649_1.fastq.gz
reads2=$wkdir/SRR8615649_2.fastq.gz
primers=$wkdir/KCNMA1_STREX_primers.fa


# Install and load software
source $HOME/.bash_profile
#conda create --name bbmap bioconda::bbmap -y
conda activate bbmap
#conda install seqtk -y


# Map reads to primers & retain only those that map

# hdist controls the number of substitutions allowed
# k sets the length of sequence that must match
# by default this also looks for the reverse compliment
# mkh = minimum kmer hits required to report the read as matching
K=21
HDIST=1
mkh=1

# map the reads to primers
bbduk.sh \
-Xmx2048M \
in=$reads1 \
in2=$reads2 \
k=$K \
hdist=$HDIST \
mkh=$mkh \
out=$wkdir/$(basename "${reads1%_*}")_without_primers_k${K}_hdist${HDIST}_mkh${mkh}.fastq \
outm=$wkdir/$(basename "${reads1%_*}")_with_primers_k${K}_hdist${HDIST}_mkh${mkh}.fastq \
stats=$wkdir/$(basename "${reads1%_*}")_stats_k${K}_hdist${HDIST}_mkh${mkh}.txt \
ref=$primers


# convert the real matches to fasta format
seqtk seq -A $wkdir/$(basename "${reads1%_*}")_with_primers_k${K}_hdist${HDIST}_mkh${mkh}.fastq > $wkdir/$(basename "${reads1%_*}")_reads_with_primers.fasta

# remove unnecessary files
rm $wkdir/$(basename "${reads1%_*}")_without_primers_k${K}_hdist${HDIST}_mkh${mkh}.fastq $wkdir/$(basename "${reads1%_*}")_with_primers_k${K}_hdist${HDIST}_mkh${mkh}.fastq



