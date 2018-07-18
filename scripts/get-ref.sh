#!/bin/bash


##
## get reference of specified type for specified genome
##


# script filename
script_name=$(basename "${BASH_SOURCE[0]}")

# check for correct number of arguments
if [ ! $# == 2 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name genome_dir ref_type \n" >&2
	exit 1
fi

# arguments
genome_root=$1
ref_type=$2


#########################


function find_file {

	local file_name=$1

	# find the shortest result
	local result=$(find -L "$genome_root" -maxdepth 2 -type f -name "$file_name" | awk '{ print length, $0 }' | sort -n | cut -d " " -f 2 | head -1)

	if [ -s "$result" ] && [ "$result" ] ; then
		echo $(readlink -f "$result")
	else
		echo -e "\n $script_name ERROR: $file_name RESULT $result DOES NOT EXIST \n" >&2
		exit 1
	fi

}

function find_dir {

	local dir_name=$1

	local result=$(find -L "$genome_root" -maxdepth 1 -type d -iname "$dir_name")

	if [ -s "$result" ] && [ "$result" ] ; then
		echo $(readlink -f "$result")
	else
		echo -e "\n $script_name ERROR: $dir_name RESULT $result DOES NOT EXIST \n" >&2
		exit 1
	fi

}

function find_basename {

	local suffix=$1

	# find the shortest result
	local result=$(find -L "$genome_root" -maxdepth 2 -type f -name "genome${suffix}" | awk '{ print length, $0 }' | sort -n | cut -d " " -f 2 | head -1)

	if [ -s "$result" ] && [ "$result" ] ; then
		result=$(readlink -f "$result")
		echo ${result/${suffix}/}
	else
		echo -e "\n $script_name ERROR: genome${suffix} RESULT $result DOES NOT EXIST \n" >&2
		exit 1
	fi

}


#########################


# genome root

if [ ! -d "$genome_root" ] ; then
	echo -e "\n $script_name ERROR: DIR $genome_root DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# file references

if [ "$ref_type" == "FASTA" ] ; then
	find_file genome.fa
fi

if [ "$ref_type" == "DICT" ] ; then
	find_file genome.dict
fi

if [ "$ref_type" == "REFFLAT" ] ; then
	find_file refFlat.txt.gz
fi

if [ "$ref_type" == "GTF" ] ; then
	find_file genes.gtf
fi

if [ "$ref_type" == "CHROMSIZES" ] ; then
	find_file chrom.sizes
fi

if [ "$ref_type" == "2BIT" ] ; then
	find_file genome.2bit
fi

if [ "$ref_type" == "FASTQSCREEN" ] ; then
	find_file fastq_screen.conf
fi

if [ "$ref_type" == "RRNAINTERVALLIST" ] ; then
	find_file rRNA.interval_list
fi


# directory references

if [ "$ref_type" == "STAR" ] ; then
	find_dir star
fi

if [ "$ref_type" == "BISMARK" ] ; then
	find_dir bismark
fi

if [ "$ref_type" == "RSEM" ] ; then
	echo "$(find_dir rsem)/ref"
fi

if [ "$ref_type" == "SALMON" ] ; then
	find_dir salmon
fi


# basename (file without suffix) references

if [ "$ref_type" == "BOWTIE1" ] ; then
	find_basename .1.ebwt
fi

if [ "$ref_type" == "BOWTIE2" ] ; then
	find_basename .1.bt2*
fi

if [ "$ref_type" == "BWA" ] ; then
	echo "$(find_basename .fa.bwt).fa"
fi


#########################



# end
