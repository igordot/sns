#!/bin/bash


# run Picard AddOrReplaceReadGroups


# script filename
script_name=$(basename "${BASH_SOURCE[0]}")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 3 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir sample_name BAM \n" >&2
	if [ $# -gt 0 ] ; then echo -e "\n ARGS: $* \n" >&2 ; fi
	exit 1
fi

# arguments
proj_dir=$(readlink -f "$1")
sample=$2
bam=$3


#########################


# check that inputs exist

if [ ! -d "$proj_dir" ] ; then
	echo -e "\n $script_name ERROR: PROJ DIR $proj_dir DOES NOT EXIST \n" >&2
	exit 1
fi

if [ ! -s "$bam" ] ; then
	echo -e "\n $script_name ERROR: BAM $bam DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# settings and files

samples_csv="${proj_dir}/samples.${segment_name}.csv"

# account for both dedup (.dd.bam) non-dedup (.bam) input BAMs
bam_base=$(basename "$bam")
bam_base=${bam_base/%.bam/}

bam_rg_dir="${proj_dir}/BAM-RG"
mkdir -p "$bam_rg_dir"
bam_rg="${bam_rg_dir}/${bam_base}.bam"

# unload all loaded modulefiles
module purge
module add default-environment


#########################


# exit if output exists already

if [ -s "$bam_rg" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	echo "${sample},${bam_rg}" >> "$samples_csv"
	exit 0
fi


#########################


# Picard

module add picard-tools/2.18.20

picard_jar="${PICARD_ROOT}/libs/picard.jar"

echo
echo " * Picard: $picard_jar "
echo " * BAM in: $bam "
echo " * BAM out: $bam_rg "
echo

bash_cmd="
java -Xms8G -Xmx8G -jar $picard_jar AddOrReplaceReadGroups \
VERBOSITY=WARNING QUIET=true VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=2500000 \
SORT_ORDER=coordinate \
CREATE_INDEX=true \
RGID=${sample} RGSM=${sample} RGLB=${sample} RGPL=ILLUMINA RGPU=LANE \
INPUT=${bam} \
OUTPUT=${bam_rg}
"
echo "CMD: $bash_cmd"
eval "$bash_cmd"

sleep 5


#########################


# check that output generated

if [ ! -s "$bam_rg" ] ; then
	echo -e "\n $script_name ERROR: BAM $bam_rg NOT GENERATED \n" >&2
	exit 1
fi


#########################


# add sample and BAM to sample sheet
echo "${sample},${bam_rg}" >> "$samples_csv"

sleep 5


#########################



# end
