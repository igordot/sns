#!/bin/bash


##
## whole genome/exome/targeted sequencing somatic single/simple/short nucleotide variant calling
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

# paths
code_dir=$(dirname $(dirname "$script_path"))
sbatch_dir="${proj_dir}/logs-sbatch"

# display settings
echo
echo " * proj_dir: $proj_dir "
echo " * sample: $sample "
echo " * code_dir: $code_dir "
echo

# specify maximum runtime for sbatch job
# SBATCHTIME=0:05:00


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

# common sbatch settings
sbatch_perf="--nodes=1 --ntasks=1 --cpus-per-task=5 --mem-per-cpu=8G"
sbatch_mail="--mail-user=${USER}@nyulangone.org --mail-type=FAIL,REQUEUE"

# Mutect2
segment_mutect2="snvs-mutect2"
bash_cmd="bash ${code_dir}/segments/${segment_mutect2}.sh $proj_dir $sample_t $bam_t $sample_n $bam_n 4"
sbatch_name="--job-name=sns.${segment_mutect2}.${sample_clean}"
sbatch_cmd="sbatch --time=48:00:00 ${sbatch_name} ${sbatch_perf} ${sbatch_mail} --export=NONE --wrap='${bash_cmd}'"
echo "CMD: $sbatch_cmd"
(eval $sbatch_cmd)

sleep 5

# Strelka
segment_strelka="snvs-strelka"
bash_cmd="bash ${code_dir}/segments/${segment_strelka}.sh $proj_dir $sample_t $bam_t $sample_n $bam_n 4"
sbatch_name="--job-name=sns.${segment_strelka}.${sample_clean}"
sbatch_cmd="sbatch --time=24:00:00 ${sbatch_name} ${sbatch_perf} ${sbatch_mail} --export=NONE --wrap='${bash_cmd}'"
echo "CMD: $sbatch_cmd"
(eval $sbatch_cmd)


#########################


date



# end
