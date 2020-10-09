#!/bin/bash


##
## chip-seq peak calling
##


# script filename
script_path="${BASH_SOURCE[0]}"
script_name=$(basename "$script_path")
route_name=${script_name/%.sh/}
echo -e "\n ========== ROUTE: $route_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 3 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir treatment_sample_name control_sample_name \n" >&2
	exit 1
fi

# standard comparison route arguments
proj_dir=$(readlink -f "$1")
sample_treatment=$2
sample_control=$3

# paths
code_dir=$(dirname $(dirname "$script_path"))
sbatch_dir="${proj_dir}/logs-sbatch"

# display settings
echo
echo " * proj_dir: $proj_dir "
echo " * treatment sample: $sample_treatment "
echo " * control sample: $sample_control "
echo " * code_dir: $code_dir "
echo

# specify maximum runtime for sbatch job
# SBATCHTIME=6:00:00


#########################


# check that inputs exist

if [ ! -d "$proj_dir" ] ; then
	echo -e "\n $script_name ERROR: dir $proj_dir does not exist \n" >&2
	exit 1
fi

# set peak type
peak_type=$(bash "${code_dir}/scripts/get-set-setting.sh" "${proj_dir}/settings.txt" EXP-PEAKS-TYPE unknown)

# set q-value cutoff to call significant regions (set to 0.05 by default)
q_value=$(bash "${code_dir}/scripts/get-set-setting.sh" "${proj_dir}/settings.txt" EXP-PEAKS-MACS-Q 0.05)

# if peak type is not known, there was a problem
if [ "$peak_type" == "unknown" ] ; then
	echo -e "\n $script_name ERROR: EXP-PEAKS-TYPE must be set in settings.txt (to 'narrow' or 'broad') \n" >&2
	exit 1
fi

# get deduplicated BAMs corresponding to sample names
segment_dedup="bam-dedup-sambamba"
bam_treatment=$(grep -s -m 1 "^${sample_treatment}," "${proj_dir}/samples.${segment_dedup}.csv" | cut -d ',' -f 2)
bam_control=$(grep -s -m 1 "^${sample_control}," "${proj_dir}/samples.${segment_dedup}.csv" | cut -d ',' -f 2)

# check if the treatment BAM exists
if [ ! -s "$bam_treatment" ] ; then
	echo -e "\n $script_name ERROR: treatment BAM $bam_treatment does not exist \n" >&2
	exit 1
fi

# if the treatment and control samples are the same, ignore the control sample
if [ "$bam_treatment" == "$bam_control" ] ; then
	bam_control=""
fi


#########################


# segments

segment_macs="peaks-macs"
bash_cmd="bash ${code_dir}/segments/${segment_macs}.sh $proj_dir $peak_type $q_value $sample_treatment $bam_treatment $bam_control"
echo "CMD: $bash_cmd"
($bash_cmd)


#########################


date



# end
