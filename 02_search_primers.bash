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
#SBATCH --time=2:00:00
#SBATCH --output=/gpfs01/home/mbzlld/code_and_scripts/slurm_out_scripts/slurm-%x-%j.out

# set variables
wkdir=/gpfs01/home/mbzlld/data/rembrandt
reads1=$wkdir/SRR8615649_1.fastq.gz
cell_lines=$wkdir/cell_lines.txt


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


# create an array of read files to loop over
accessions=($(tail -n +2 $cell_lines | cut -f2))
reads=( "${accessions[@]/#/$wkdir\/}" )
reads=( "${reads[@]/%/_1.fastq.gz}" )


# create an array containing the primer files to loop over
primers=($wkdir/KCNMA1_STREX_primers.fa $wkdir/KCNMA1_primers.fa)



for file in "${reads[@]}" ; do
	for primer in "${primers[@]}" ; do
		# map the reads to primers
		bbduk.sh \
		-Xmx2048M \
		in=$file \
		k=$K \
		hdist=$HDIST \
		mkh=$mkh \
		out=$wkdir/$(basename "${file%.*.*}")_without_primers_$(basename "${primer%_*}").fastq \
		outm=$wkdir/$(basename "${file%.*.*}")_with_primers_$(basename "${primer%_*}").fastq \
		stats=$wkdir/$(basename "${file%.*.*}")_$(basename "${primer%_*}")_stats.txt \
		ref=$primer
		
		# convert the real matches to fasta format
		#seqtk seq -A $wkdir/$(basename "${file%.*.**}")_with_primers_$(basename "${primer%_*}").fastq > $wkdir/$(basename "${file%.*.*}")_reads_with_primers_$(basename "${primer%_*}").fasta
		
		# remove unnecessary files
		rm $wkdir/$(basename "${file%.*.*}")_without_primers_$(basename "${primer%_*}").fastq
		rm $wkdir/$(basename "${file%.*.*}")_with_primers_$(basename "${primer%_*}").fastq
		
	done
done

