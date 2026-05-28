#!/bin/bash


##
## Identify species corresponding to the sequencing reads
##


# specify maximum runtime for sbatch job
# SBATCHTIME=6:00:00

# standard route header (validate args, print settings, prepare environment)
code_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
source "${code_dir}/scripts/route-header.sh" "$@"


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

# check that sbatch can be run
if ! sbatch --version >/dev/null 2>&1; then
	echo -e "\n $script_name ERROR: sbatch cannot be executed \n" >&2
	exit 1
fi

# species using Centrifuge (a scan of the NCBI nt database)
# separate job due to memory requirements
# out-of-memory error with 150G (3/2024)
segment_species="species-centrifuge"
bash_cmd="bash ${code_dir}/segments/${segment_species}.sh $proj_dir $sample 4 $fastq_R1"
sbatch_perf="--nodes=1 --ntasks=1 --cpus-per-task=5 --mem=200G"
sbatch_mail="--mail-user=${USER}@nyulangone.org --mail-type=FAIL,REQUEUE"
sbatch_name="--job-name=sns.${segment_species}.${sample}"
sbatch_cmd="sbatch --time=4:00:00 ${sbatch_name} ${sbatch_perf} ${sbatch_mail} --export=NONE --wrap='${bash_cmd}'"
echo "CMD: $sbatch_cmd"
(eval $sbatch_cmd)

# species using FastQ Screen (a quick scan of a small set of pre-defined common species)
segment_species="species-fastqscreen"
bash_cmd="bash ${code_dir}/segments/${segment_species}.sh $proj_dir $sample $threads $fastq_R1"
($bash_cmd)


#########################


date



# end
