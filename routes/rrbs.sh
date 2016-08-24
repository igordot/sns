#!/bin/bash


##
## RRBS using Bismark
##


# script filename
script_name=$(basename "${BASH_SOURCE[0]}")
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
code_dir=$(dirname "$(dirname "${BASH_SOURCE[0]}")")
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
fastq_R1=$(grep -s -m 1 "^${sample}," "${proj_dir}/samples.fastq-clean.csv" | cut -d ',' -f 2)
fastq_R2=$(grep -s -m 1 "^${sample}," "${proj_dir}/samples.fastq-clean.csv" | cut -d ',' -f 3)
if [ -z "$fastq_R1" ] ; then
	bash_cmd="bash ${code_dir}/segments/fastq-fastq-clean.sh $proj_dir $sample"
	($bash_cmd)
	fastq_R1=$(grep -m 1 "^${sample}," "${proj_dir}/samples.fastq-clean.csv" | cut -d ',' -f 2)
	fastq_R2=$(grep -m 1 "^${sample}," "${proj_dir}/samples.fastq-clean.csv" | cut -d ',' -f 3)
fi

# trim FASTQs with Trim Galore
fastq_R1_trimmed=$(grep -s -m 1 "^${sample}," "${proj_dir}/samples.fastq-trim.csv" | cut -d ',' -f 2)
fastq_R2_trimmed=$(grep -s -m 1 "^${sample}," "${proj_dir}/samples.fastq-trim.csv" | cut -d ',' -f 3)
if [ -z "$fastq_R1_trimmed" ] ; then
	bash_cmd="bash ${code_dir}/segments/fastq-fastq-trim-trimgalore.sh $proj_dir $sample rrbs $fastq_R1 $fastq_R2"
	($bash_cmd)
	fastq_R1_trimmed=$(grep -m 1 "^${sample}," "${proj_dir}/samples.fastq-trim.csv" | cut -d ',' -f 2)
	fastq_R2_trimmed=$(grep -m 1 "^${sample}," "${proj_dir}/samples.fastq-trim.csv" | cut -d ',' -f 3)
fi

# run Bismark alignment
bam_bismark=$(grep -s -m 1 "^${sample}," "${proj_dir}/samples.bam-bismark.csv" | cut -d ',' -f 2)
if [ -z "$bam_bismark" ] ; then
	bash_cmd="bash ${code_dir}/segments/fastq-bam-bismark.sh $proj_dir $sample $threads $fastq_R1_trimmed $fastq_R2_trimmed"
	($bash_cmd)
	bam_bismark=$(grep -m 1 "^${sample}," "${proj_dir}/samples.bam-bismark.csv" | cut -d ',' -f 2)
fi

# run Bismark methylation extractor
if [ -n "$fastq_R2" ] ; then
	#
	echo "pe"
else
	bash_cmd="bash ${code_dir}/segments/bam-meth-bismark.sh $proj_dir $sample $threads $bam_bismark se"
	($bash_cmd)
	bash_cmd="bash ${code_dir}/segments/bam-meth-bismark.sh $proj_dir $sample $threads $bam_bismark se-ignore-r1-3"
	($bash_cmd)
fi


#########################


# delete empty qsub .po files
rm -f ${qsub_dir}/sns.*.po*


#########################


date



# end
