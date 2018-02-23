#!/bin/bash


# remove duplicates using sambamba


# script filename
script_name=$(basename "${BASH_SOURCE[0]}")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 4 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir sample_name num_threads BAM \n" >&2
	if [ $# -gt 0 ] ; then echo -e "\n ARGS: $* \n" >&2 ; fi
	exit 1
fi

# arguments
proj_dir=$1
sample=$2
threads=$3
bam=$4


#########################


# check that inputs exist

if [ ! -d "$proj_dir" ] ; then
	echo -e "\n $script_name ERROR: DIR $proj_dir DOES NOT EXIST \n" >&2
	exit 1
fi

if [ ! -s "$bam" ] ; then
	echo -e "\n $script_name ERROR: BAM $bam DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# settings and files

summary_dir="${proj_dir}/summary"
mkdir -p "$summary_dir"
summary_csv="${summary_dir}/${sample}.${segment_name}.csv"

samples_csv="${proj_dir}/samples.${segment_name}.csv"

bam_dd_dir="${proj_dir}/BAM-DD"
mkdir -p "$bam_dd_dir"
bam_dd="${bam_dd_dir}/${sample}.dd.bam"
bai_dd="${bam_dd}.bai"

dedup_logs_dir="${proj_dir}/logs-${segment_name}"
mkdir -p "$dedup_logs_dir"
bam_dd_log="${dedup_logs_dir}/${sample}.log.txt"
bam_dd_flagstat="${dedup_logs_dir}/${sample}.flagstat.txt"

# unload all loaded modulefiles
module purge
module load local


#########################


# exit if output exists already

if [ -s "$bam_dd" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	echo -e "\n $script_name ADD $sample TO $samples_csv \n" >&2
	echo "${sample},${bam_dd}" >> "$samples_csv"
	exit 1
fi


#########################


# sambamba markdup

sambamba_bin="/ifs/home/id460/software/sambamba/sambamba_v0.6.7"

echo " * sambamba: $(readlink -f $(which $sambamba_bin)) "
echo " * sambamba version: $($sambamba_bin 2>&1 | head -1) "
echo " * BAM IN: $bam "
echo " * BAM OUT: $bam_dd "

bash_cmd="
$sambamba_bin markdup \
--remove-duplicates \
--nthreads $threads \
--hash-table-size 525000 \
--overflow-list-size 525000 \
$bam \
$bam_dd \
2> $bam_dd_log
"
echo "CMD: $bash_cmd"
eval "$bash_cmd"

sleep 30


#########################


# check that sambamba markdup output generated

if [ ! -s "$bam_dd" ] ; then
	echo -e "\n $script_name ERROR: BAM $bam_dd NOT GENERATED \n" >&2
	exit 1
fi

if [ ! -s "$bai_dd" ] ; then
	echo -e "\n $script_name ERROR: BAI $bai_dd NOT GENERATED \n" >&2
	# delete BAM since something went wrong and it might be corrupted
	rm -fv "$bam_dd"
	exit 1
fi


#########################


# run flagstat

bash_cmd="$sambamba_bin flagstat $bam_dd > $bam_dd_flagstat"
echo "CMD: $bash_cmd"
eval "$bash_cmd"

sleep 30


#########################


# check that flagstat output generated

if [ ! -s "$bam_dd_flagstat" ] ; then
	echo -e "\n $script_name ERROR: FLAGSTAT $bam_dd_flagstat NOT GENERATED \n" >&2
	# delete BAM and BAI since something went wrong and they might be corrupted
	rm -fv "$bam_dd"
	rm -fv "$bai_dd"
	exit 1
fi


#########################


# generate alignment summary

reads_duplicates=$(cat "$bam_dd_log" | grep -m 1 "found.*duplicates" | tr -d -c 0-9)

# "mapped" and "total" lines include secondary alignments
reads_dedup_all=$(cat "$bam_dd_flagstat" | grep -m 1 "mapped (" | cut -d ' ' -f 1)
reads_dedup_sec=$(cat "$bam_dd_flagstat" | grep -m 1 "secondary" | cut -d ' ' -f 1)
reads_dedup=$(echo "${reads_dedup_all} - ${reads_dedup_sec}" | bc)

reads_mapped=$(echo "${reads_dedup} + ${reads_duplicates}" | bc)

reads_duplicates_pct=$(echo "(${reads_duplicates} / ${reads_mapped}) * 100" | bc -l | cut -c 1-4)
reads_duplicates_pct="${reads_duplicates_pct}%"

# header for summary file
echo "#SAMPLE,MAPPED READS,DEDUPLICATED READS,DUPLICATES %" > "$summary_csv"

# summarize log file
echo "${sample},${reads_mapped},${reads_dedup},${reads_duplicates_pct}" >> "$summary_csv"

sleep 30

# combine all sample summaries
cat ${summary_dir}/*.${segment_name}.csv | LC_ALL=C sort -t ',' -k1,1 | uniq > "${proj_dir}/summary.${segment_name}.csv"


#########################


# add sample and BAM to sample sheet
echo "${sample},${bam_dd}" >> "$samples_csv"

sleep 30


#########################



# end
