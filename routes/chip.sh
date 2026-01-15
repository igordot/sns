#!/bin/bash


##
## ChIP-seq using Bowtie 2
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

# paths
code_dir=$(dirname $(dirname "$script_path"))

# reserve a thread for overhead
threads=$SLURM_CPUS_PER_TASK
threads=$(( threads - 1 ))

# show settings
echo
echo " * proj_dir: $proj_dir "
echo " * sample: $sample "
echo " * code_dir: $code_dir "
echo " * slurm threads: $SLURM_CPUS_PER_TASK "
echo " * command threads: $threads "
echo " * slurm nodename: $SLURMD_NODENAME "
echo " * hostname: $(hostname) "
echo " * time: $(date "+%Y-%m-%d %H:%M") "
echo

# specify maximum runtime for sbatch job
# SBATCHTIME=48:00:00


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

# run alignment
segment_align="align-bowtie2-chip"
bam_bt2=$(grep -s -m 1 "^${sample}," "${proj_dir}/samples.${segment_align}.csv" | cut -d ',' -f 2)
if [ -z "$bam_bt2" ] ; then
	bash_cmd="bash ${code_dir}/segments/${segment_align}.sh $proj_dir $sample $threads $fastq_R1 $fastq_R2"
	($bash_cmd)
	bam_bt2=$(grep -m 1 "^${sample}," "${proj_dir}/samples.${segment_align}.csv" | cut -d ',' -f 2)
fi

# if BAM is not set, there was a problem
if [ -z "$bam_bt2" ] ; then
	echo -e "\n $script_name ERROR: SEGMENT $segment_align DID NOT FINISH \n" >&2
	exit 1
fi

# remove duplicates
segment_dedup="bam-dedup-sambamba"
bam_dd=$(grep -s -m 1 "^${sample}," "${proj_dir}/samples.${segment_dedup}.csv" | cut -d ',' -f 2)
if [ -z "$bam_dd" ] ; then
	bash_cmd="bash ${code_dir}/segments/${segment_dedup}.sh $proj_dir $sample $threads $bam_bt2"
	($bash_cmd)
	bam_dd=$(grep -m 1 "^${sample}," "${proj_dir}/samples.${segment_dedup}.csv" | cut -d ',' -f 2)
fi

# if BAM is not set, there was a problem
if [ -z "$bam_dd" ] ; then
	echo -e "\n $script_name ERROR: SEGMENT $segment_dedup DID NOT FINISH \n" >&2
	exit 1
fi

# generate BigWig (deeptools)
segment_bw_deeptools="bigwig-deeptools"
bash_cmd="bash ${code_dir}/segments/${segment_bw_deeptools}.sh $proj_dir $sample 4 $bam_dd"
sbatch_perf="--nodes=1 --ntasks=1 --cpus-per-task=5 --mem=25G"
sbatch_mail="--mail-user=${USER}@nyulangone.org --mail-type=FAIL,REQUEUE"
sbatch_name="--job-name=sns.${segment_bw_deeptools}.${sample}"
sbatch_cmd="sbatch --time=8:00:00 ${sbatch_name} ${sbatch_perf} ${sbatch_mail} --export=NONE --wrap='${bash_cmd}'"
echo "CMD: $sbatch_cmd"
(eval $sbatch_cmd)

# set peak type for peak calling (should be manually adjusted to the proper value)
bash "${code_dir}/scripts/get-set-setting.sh" "${proj_dir}/settings.txt" EXP-PEAKS-TYPE unknown


#########################


# combine summary from each step

sleep 5

summary_csv="${proj_dir}/summary-combined.${route_name}.csv"

bash_cmd="
bash ${code_dir}/scripts/join-many.sh , X \
${proj_dir}/summary.${segment_fastq_clean}.csv \
${proj_dir}/summary.${segment_align}.csv \
${proj_dir}/summary.${segment_dedup}.csv \
> $summary_csv
"
(eval $bash_cmd)


#########################


# generate pairs sample sheet template

samples_pairs_csv="${proj_dir}/samples.pairs.csv"

if [ ! -s "$samples_pairs_csv" ] ; then
	echo "#SAMPLE-TREATMENT,#SAMPLE-CONTROL" > $samples_pairs_csv
	sed 's/\,.*/,NA/g' ${proj_dir}/samples.fastq-raw.csv | LC_ALL=C sort -u >> $samples_pairs_csv
fi


#########################


date



# end
