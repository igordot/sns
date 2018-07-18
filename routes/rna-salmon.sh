#!/bin/bash


##
## RNA-seq using Salmon transcript quantification
##


# script filename
script_path="${BASH_SOURCE[0]}"
script_name=$(basename "$script_path")
route_name=${script_name/%.sh/}
echo -e "\n ========== ROUTE: $route_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 2 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir sample_name \n" >&2
	exit 1
fi

# standard route arguments
proj_dir=$(readlink -f "$1")
sample=$2

# additional settings
threads=$NSLOTS
code_dir=$(dirname $(dirname "$script_path"))
qsub_dir="${proj_dir}/logs-qsub"

# display settings
echo " * proj_dir: $proj_dir "
echo " * sample: $sample "
echo " * code_dir: $code_dir "
echo " * qsub_dir: $qsub_dir "
echo " * threads: $threads "


#########################


# delete empty qsub .po files
rm -f ${qsub_dir}/sns.*.po*


#########################


# segments

# rename and/or merge raw input FASTQs
segment_fastq_clean="fastq-clean"
fastq_R1=$(grep -s -m 1 "^${sample}," "${proj_dir}/samples.${segment_fastq_clean}.csv" | cut -d ',' -f 2)
fastq_R2=$(grep -s -m 1 "^${sample}," "${proj_dir}/samples.${segment_fastq_clean}.csv" | cut -d ',' -f 3)
if [ -z "$fastq_R1" ] ; then
	bash_cmd="bash ${code_dir}/segments/${segment_fastq_clean}.sh $proj_dir $sample"
	($bash_cmd)
	fastq_R1=$(grep -m 1 "^${sample}," "${proj_dir}/samples.${segment_fastq_clean}.csv" | cut -d ',' -f 2)
	fastq_R2=$(grep -m 1 "^${sample}," "${proj_dir}/samples.${segment_fastq_clean}.csv" | cut -d ',' -f 3)
fi

# if FASTQ is not set, there was a problem
if [ -z "$fastq_R1" ] ; then
	echo -e "\n $script_name ERROR: SEGMENT $segment_fastq_clean DID NOT FINISH \n" >&2
	exit 1
fi

# run FastQC (separately for paired-end reads)
segment_qc_fastqc="qc-fastqc"
bash_cmd="bash ${code_dir}/segments/${segment_qc_fastqc}.sh $proj_dir $sample $threads $fastq_R1"
($bash_cmd)
if [ -n "$fastq_R2" ] ; then
	bash_cmd="bash ${code_dir}/segments/${segment_qc_fastqc}.sh $proj_dir $sample $threads $fastq_R2"
	($bash_cmd)
fi

# fastq_screen
# bash_cmd="bash ${code_dir}/segments/qc-fastqscreen.sh $proj_dir $sample $fastq_R1"
# ($bash_cmd)

# Salmon
segment_quant="quant-salmon"
bash_cmd="bash ${code_dir}/segments/${segment_quant}.sh $proj_dir $sample $threads $fastq_R1 $fastq_R2"
($bash_cmd)


#########################


# combine summary from each step

sleep 30

summary_csv="${proj_dir}/summary-combined.${route_name}.csv"

bash_cmd="
bash ${code_dir}/scripts/join-many.sh , X \
${proj_dir}/summary.${segment_fastq_clean}.csv \
${proj_dir}/summary.${segment_quant}.csv \
> $summary_csv
"
(eval $bash_cmd)


#########################


# delete empty qsub .po files
rm -f ${qsub_dir}/sns.*.po*


#########################


date



# end
