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
#SBATCH --time=10:00:00
#SBATCH --output=/gpfs01/home/mbzlld/code_and_scripts/slurm_out_scripts/slurm-%x-%j.out

# set variables
wkdir=/gpfs01/home/mbzlld/data/rembrandt

#cell_lines=$wkdir/cell_lines.txt
#output_file=results.txt

cell_lines=$wkdir/SF188_accessions.txt
output_file=SF188_results.txt

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


# loop over the read files for each cell line 
for file in "${reads[@]}" ; do
	# loop over each primer
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

# deactivate software
conda deactivate

#####################################
# Generate a summary of the results #
#####################################


############## FOR THE STREX OUTPUT FILES ###############

# generate an array of output files to loop over
outputs_STREX=( "${accessions[@]/#/$wkdir\/}" )
outputs_STREX=( "${outputs_STREX[@]/%/_1_KCNMA1_STREX_stats.txt}" )


# loop over output stats files and add to results
: > $wkdir/$output_file # create empty file to write to (empty the file if it already exists)
i=1 # create for loop numbering
for file in ${outputs_STREX[@]} ; do
	
	# start results file with list of accessions from output stats files
	accession_no=$( sed -n '1p' $file | cut -f2 )
	accession_no=$( basename ${accession_no%_*} )
	echo $accession_no >> $wkdir/$output_file
	
	# add new column of total reads info to results file
	sed -i "${i}s/$/\t$( sed -n '2p' $file | cut -f2 )/" $wkdir/$output_file
	
	# add new column of STREX splice variant reads to output
	sed -i "${i}s/$/\t$( sed -n '3p' $file | cut -f2 )/" $wkdir/$output_file
	
	((i++))
done


############## FOR THE KCNMA1 OUTPUT FILES ###############

# generate an array of output files to loop over
outputs_KCNMA1=( "${accessions[@]/#/$wkdir\/}" )
outputs_KCNMA1=( "${outputs_KCNMA1[@]/%/_1_KCNMA1_stats.txt}" )


# loop over output stats files and add to results
i=1 # create for loop numbering
for file in ${outputs_KCNMA1[@]} ; do

        # add new column of KCNM1 reads to output
        sed -i "${i}s/$/\t$( sed -n '3p' $file | cut -f2 )/" $wkdir/$output_file
	
        ((i++))
done

######################

# Finally add a header line
sed -i '1s/^/accession\ttotal_reads\tKCNMA1_STREX_reads\tKCNMA1_reads\n/' $wkdir/$output_file

# and add percentage calculations
#awk 'BEGIN {OFS="\t"} NR==1 {print $0, "percent_STREX"} NR>1 {printf "%s\t%.10f\n", $0, ($3 / $2) * 100}' $wkdir/results.txt > $wkdir/tmp && mv $wkdir/tmp $wkdir/results.txt
awk 'BEGIN {OFS="\t"} NR==1 {print $0, "percent_STREX"} 
NR>1 {
  if ($2 == 0) {
    print $0, "NA"
  } else {
    printf "%s\t%.10f\n", $0, ($3 / $2) * 100
  }
} ' $wkdir/$output_file > $wkdir/tmp && mv $wkdir/tmp $wkdir/$output_file



#awk 'BEGIN {OFS="\t"} NR==1 {print $0, "percent_KCNMA1"} NR>1 {printf "%s\t%.10f\n", $0, ($4 / $2) * 100}' $wkdir/results.txt > $wkdir/tmp && mv $wkdir/tmp $wkdir/results.txt
awk 'BEGIN {OFS="\t"} NR==1 {print $0, "percent_KCNMA1"} 
NR>1 {
  if ($2 == 0) {
	  print $0, "NA"
  } else {
  printf "%s\t%.10f\n", $0, ($4 / $2) * 100}}' $wkdir/$output_file > $wkdir/tmp && mv $wkdir/tmp $wkdir/$output_file

#awk 'BEGIN {OFS="\t"} NR==1 {print $0, "proportion_KCNMA1_STREX"} NR>1 {printf "%s\t%.10f\n", $0, ($5 / $6) * 100}' $wkdir/results.txt > $wkdir/tmp && mv $wkdir/tmp $wkdir/results.txt
awk 'BEGIN {OFS="\t"} NR==1 {print $0, "proportion_KCNMA1_STREX"} NR>1 { if ($6 == 0) { print $0, "NA" } else { printf "%s\t%.10f\n", $0, ($5 / $6) * 100}}' $wkdir/$output_file > $wkdir/tmp && mv $wkdir/tmp $wkdir/

# add cell line column to the results file from the cell_lines file
awk -F'\t' '
    FNR==NR { map[$2]=$1; next }           # Read file2.tsv into a map using col1 as key and col2 as value
    FNR==1 { print "cell_line", $0; next } # Add header to output
    { print map[$1], $0 }                  # Prepend matched value from file2.tsv using col2 in file1.tsv as key
' OFS='\t' $cell_lines $wkdir/$output_file > $wkdir/tmp && mv $wkdir/tmp $wkdir/$output_file

######################


