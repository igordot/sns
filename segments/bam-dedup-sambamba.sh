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
proj_dir=$(readlink -f "$1")
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
bam_dd_temp="${bam_dd_dir}/${sample}.dd.tmp.bam"
bam_dd_final="${bam_dd_dir}/${sample}.dd.bam"

dedup_logs_dir="${proj_dir}/logs-${segment_name}"
mkdir -p "$dedup_logs_dir"
bam_dd_log="${dedup_logs_dir}/${sample}.log.txt"
bam_dd_flagstat="${dedup_logs_dir}/${sample}.flagstat.txt"

# unload all loaded modulefiles
module purge
module add default-environment


#########################


# check for output

# skip if final BAM and BAI exist
if [ -s "$bam_dd_final" ] && [ -s "${bam_dd_final}.bai" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	echo -e "\n $script_name ADD $sample TO $samples_csv \n" >&2
	echo "${sample},${bam_dd_final}" >> "$samples_csv"
	exit 0
fi

# delete temp BAM (likely incomplete since the final BAM was not generated)
if [ -s "$bam_dd_temp" ] ; then
	echo -e "\n $script_name WARNING: POTENTIALLY CORRUPT BAM $bam_dd_temp EXISTS \n" >&2
	rm -fv "$bam_dd_temp"
	rm -fv "${bam_dd_temp}.bai"
fi


#########################


# sambamba markdup

module add sambamba/0.6.8

sambamba_bin="sambamba-0.6.8"

echo
echo " * sambamba: $(readlink -f $(which $sambamba_bin)) "
echo " * sambamba version: $($sambamba_bin 2>&1 | grep -m 1 'sambamba') "
echo " * BAM in: $bam "
echo " * BAM out temp: $bam_dd_temp "
echo " * BAM out final: $bam_dd_final "
echo

bash_cmd="
$sambamba_bin markdup \
--remove-duplicates \
--nthreads $threads \
--hash-table-size 525000 \
--overflow-list-size 525000 \
$bam \
$bam_dd_temp \
2> $bam_dd_log
"
echo "CMD: $bash_cmd"
eval "$bash_cmd"

sleep 5


#########################


# check that sambamba markdup output generated

if [ ! -s "$bam_dd_temp" ] ; then
	echo -e "\n $script_name ERROR: BAM $bam_dd_temp NOT GENERATED \n" >&2
	exit 1
fi

if [ ! -s "${bam_dd_temp}.bai" ] ; then
	echo -e "\n $script_name ERROR: BAI ${bam_dd_temp}.bai NOT GENERATED \n" >&2
	# delete BAM since something went wrong and it might be corrupted
	rm -fv "$bam_dd_temp"
	exit 1
fi


#########################


# run flagstat

bash_cmd="$sambamba_bin flagstat $bam_dd_temp > $bam_dd_flagstat"
echo "CMD: $bash_cmd"
eval "$bash_cmd"

sleep 5


#########################


# check that flagstat output generated

if [ ! -s "$bam_dd_flagstat" ] ; then
	echo -e "\n $script_name ERROR: FLAGSTAT $bam_dd_flagstat NOT GENERATED \n" >&2
	# delete BAM and BAI since something went wrong and they might be corrupted
	rm -fv "$bam_dd_temp"
	rm -fv "${bam_dd_temp}.bai"
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

sleep 5

# combine all sample summaries
cat ${summary_dir}/*.${segment_name}.csv | LC_ALL=C sort -t ',' -k1,1 | uniq > "${proj_dir}/summary.${segment_name}.csv"


#########################


# finalize

# move BAM and BAI from temp to final
mv -v "$bam_dd_temp" "${bam_dd_final}"
mv -v "${bam_dd_temp}.bai" "${bam_dd_final}.bai"

# add sample and BAM to sample sheet
echo "${sample},${bam_dd_final}" >> "$samples_csv"

sleep 5


#########################



# end
