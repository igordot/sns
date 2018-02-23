#!/bin/bash


# run FastQC


# script filename
script_path="${BASH_SOURCE[0]}"
script_name=$(basename "$script_path")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 4 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir sample_name num_threads FASTQ \n" >&2
	if [ $# -gt 0 ] ; then echo -e "\n ARGS: $* \n" >&2 ; fi
	exit 1
fi

# arguments
proj_dir=$1
sample=$2
threads=$3
fastq=$4


#########################


# settings and files

out_dir="${proj_dir}/QC-FastQC"
mkdir -p "$out_dir"
fastqc_basename=$(basename "$fastq")
fastqc_basename=${fastqc_basename/%.fastq.gz/_fastqc}
report_html="${out_dir}/${fastqc_basename}.html"
report_zip="${out_dir}/${fastqc_basename}.zip"

# unload all loaded modulefiles
module purge
module load local


#########################


# exit if output exits already

if [ -s "$report_html" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 1
fi


#########################


# check that inputs exist

if [ ! -d "$proj_dir" ] || [ ! "$proj_dir" ] ; then
	echo -e "\n $script_name ERROR: DIR $proj_dir DOES NOT EXIST \n" >&2
	exit 1
fi

if [ ! -s "$fastq" ] || [ ! "$fastq" ] ; then
	echo -e "\n $script_name ERROR: FASTQ $fastq DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# run FastQC

module load fastqc/0.11.7

echo
echo " * FastQC path: $(readlink -f $(which fastqc)) "
echo " * FastQC version: $(fastqc --version) "
echo " * FASTQ: $fastq "
echo " * report html: $report_html "
echo " * report zip: $report_zip "
echo

# --quiet      supress all progress messages on stdout and only report errors
# --noextract  do not uncompress the output file after creating it
# --nogroup    disable grouping of bases for reads >50bp

fastqc_cmd="fastqc \
--quiet --nogroup --noextract \
--threads $threads \
--outdir $out_dir \
$fastq
"
echo "CMD: $fastqc_cmd"
$fastqc_cmd


#########################


# check that output generated

if [ ! -s "$report_html" ] ; then
	echo -e "\n $script_name ERROR: REPORT $report_html NOT GENERATED \n" >&2
	exit 1
fi

if [ ! -s "$report_zip" ] ; then
	echo -e "\n $script_name ERROR: REPORT $report_zip NOT GENERATED \n" >&2
	exit 1
fi


#########################



# end
