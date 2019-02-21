#!/bin/bash


# labeling of reads and quantification of species with Centrifuge


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


# check that inputs exist

if [ ! -d "$proj_dir" ] ; then
	echo -e "\n $script_name ERROR: DIR $proj_dir DOES NOT EXIST \n" >&2
	exit 1
fi

if [ ! -s "$fastq" ] ; then
	echo -e "\n $script_name ERROR: FASTQ $fastq DOES NOT EXIST \n" >&2
	exit 1
fi

centrifuge_index="/gpfs/data/igorlab/ref/Centrifuge/nt_2018_3_3/nt"

if [ ! -s "${centrifuge_index}.1.cf" ] ; then
	echo -e "\n $script_name ERROR: REF $centrifuge_index DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# settings and files

# summary_dir="${proj_dir}/summary"
# mkdir -p "$summary_dir"
# summary_csv="${summary_dir}/${sample}.${segment_name}.csv"

centrifuge_dir="${proj_dir}/species-centrifuge"
mkdir -p "$centrifuge_dir"
report_csv="${centrifuge_dir}/${sample}.csv"

logs_dir="${proj_dir}/logs-${segment_name}"
mkdir -p "$logs_dir"
results_txt="${logs_dir}/${sample}.results.txt"
report_txt="${logs_dir}/${sample}.report.txt"

# unload all loaded modulefiles
module purge
module add default-environment


#########################


# exit if output exits already

if [ -s "$report_csv" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 0
fi


#########################


# Centrifuge

centrifuge_bin="/gpfs/data/igorlab/software/Centrifuge/centrifuge-1.0.4-beta/bin/centrifuge"

echo
echo " * Centrifuge: $(readlink -f $(which $centrifuge_bin)) "
echo " * Centrifuge version: $($centrifuge_bin --version 2>&1 | head -1) "
echo " * FASTQ: $fastq "
echo " * results: $results_txt "
echo " * report: $report_txt "
echo

# using 1M reads
# not using metrics file because it is not very helpful

bash_cmd="
zcat $fastq \
| head -4000000 \
| $centrifuge_bin \
--threads $threads \
-x $centrifuge_index \
-q -U - \
--report-file $report_txt \
-S $results_txt
"
echo "CMD: $bash_cmd"
eval "$bash_cmd"

sleep 5


#########################


# check that output generated

if [ ! -s "$report_txt" ] ; then
	echo -e "\n $script_name ERROR: REPORT $report_txt NOT GENERATED \n" >&2
	exit 1
fi


#########################


# clean up results report

# col 1: name of a genome
# col 2: taxonomic ID (e.g., 36870)
# col 3: taxonomic rank (e.g., leaf)
# col 4: length of the genome sequence (e.g., 703004)
# col 5: number of reads classified to this genomic sequence including multi-classified reads (e.g., 5981)
# col 6: number of reads uniquely classified to this genomic sequence (e.g., 5964)
# col 7: proportion of this genome normalized by its genomic length (e.g., 0.0152317)

# genomeSize and abundance can be incorrect (https://github.com/infphilo/centrifuge/issues/13)

# header
cat "$report_txt" | grep -m 1 "numUniqueReads" | cut -f 1,2,3,5,6 | tr '\t' ',' > $report_csv

# most abundant species by unique reads
cat "$report_txt" \
| grep -v "numUniqueReads" \
| sort -t $'\t' -k6,6nr \
| head -20 \
| cut -f 1,2,3,5,6 \
| tr ',' ' ' \
| tr '\t' ',' \
>> "$report_csv"

# delete classification results file (contains info for each input read)
rm -fv "$results_txt"


#########################



# end
