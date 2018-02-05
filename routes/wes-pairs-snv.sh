#!/bin/bash


##
## whole genome/exome/targeted sequencing somatic single/short nucleotide variant calling
##


# script filename
script_path="${BASH_SOURCE[0]}"
script_name=$(basename "$script_path")
route_name=${script_name/%.sh/}
echo -e "\n ========== ROUTE: $route_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 3 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir tumor_sample_name normal_sample_name \n" >&2
	exit 1
fi

# standard comparison route arguments
proj_dir=$(readlink -f "$1")
sample_t=$2
sample_n=$3
sample_clean="${sample_t}-${sample_n}"

# additional settings
code_dir=$(dirname $(dirname "$script_path"))
qsub_dir="${proj_dir}/logs-qsub"

# display settings
echo
echo " * proj_dir: $proj_dir "
echo " * sample T: $sample_t "
echo " * sample N: $sample_n "
echo " * code_dir: $code_dir "
echo " * qsub_dir: $qsub_dir "
echo


#########################


# delete empty qsub .po files
rm -f ${qsub_dir}/sns.*.po*


#########################


# check that inputs exist

if [ ! -d "$proj_dir" ] ; then
	echo -e "\n $script_name ERROR: DIR $proj_dir DOES NOT EXIST \n" >&2
	exit 1
fi

# get GATK-processed BAMs corresponding to sample names
segment_gatk="bam-ra-rc-gatk"
bam_t=$(grep -s -m 1 "^${sample_t}," "${proj_dir}/samples.${segment_gatk}.csv" | cut -d ',' -f 2)
bam_n=$(grep -s -m 1 "^${sample_n}," "${proj_dir}/samples.${segment_gatk}.csv" | cut -d ',' -f 2)

if [ ! -s "$bam_t" ] ; then
	echo -e "\n $script_name ERROR: BAM $bam_t DOES NOT EXIST \n" >&2
	exit 1
fi

if [ ! -s "$bam_n" ] ; then
	echo -e "\n $script_name ERROR: BAM $bam_n DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# segments

qsub_settings="-q all.q -M ${USER}@nyumc.org -m a -j y -cwd -b y -hard -l mem_free=150G"

# MuTect2
segment_mutect2="snvs-mutect2"
bash_cmd="bash ${code_dir}/segments/${segment_mutect2}.sh $proj_dir $sample_t $bam_t $sample_n $bam_n"
qsub_cmd="qsub -N sns.${segment_mutect2}.${sample_clean} ${qsub_settings} ${bash_cmd}"
echo "CMD: $qsub_cmd"
($qsub_cmd)

sleep 30

# Strelka
segment_strelka="snvs-strelka"
bash_cmd="bash ${code_dir}/segments/${segment_strelka}.sh $proj_dir $sample_t $bam_t $sample_n $bam_n 8"
qsub_cmd="qsub -N sns.${segment_strelka}.${sample_clean} ${qsub_settings} -pe threaded 8 ${bash_cmd}"
echo "CMD: $qsub_cmd"
($qsub_cmd)


#########################


# delete empty qsub .po files
rm -f ${qsub_dir}/sns.*.po*


#########################


date



# end
