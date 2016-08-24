#!/bin/bash


# Bismark deduplicate_bismark


# script filename
script_name=$(basename "${BASH_SOURCE[0]}")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 4 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name [project dir] [sample] [Bismark BAM] [analysis type] \n" >&2
	exit 1
fi

# arguments
proj_dir=$1
sample=$2
bismark_bam=$3
analysis_type=$4


#########################


# check that inputs exist

if [ ! -d "$proj_dir" ] ; then
	echo -e "\n $script_name ERROR: DIR $proj_dir DOES NOT EXIST \n" >&2
	exit 1
fi

if [ ! -s "$bismark_bam" ] ; then
	echo -e "\n $script_name ERROR: BAM $bismark_bam DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# settings and files

summary_dir="${proj_dir}/summary"
mkdir -p "$summary_dir"
summary_csv="${summary_dir}/${sample}.${segment_name}.csv"

samples_csv="${proj_dir}/samples.${segment_name}.csv"

bismark_bam_dd_dir="${proj_dir}/BAM-DD-Bismark"
mkdir -p "$bismark_bam_dd_dir"
bismark_bam_dd_final="${bismark_bam_dd_dir}/${sample}.bam"

bismark_dd_logs_dir="${proj_dir}/logs-${segment_name}"
mkdir -p "$bismark_dd_logs_dir"
bismark_dd_report_final="${bismark_dd_logs_dir}/${sample}.report.txt"

bismark_report_dir="${proj_dir}/Bismark-report"
mkdir -p "$bismark_report_dir"
bismark_report_deduplication="${bismark_report_dir}/${sample}_bismark_bt2_pe.deduplication_report.txt"

bismark_bam_dd_original="${bismark_bam/%.bam/.deduplicated.bam}"
bismark_dd_report_original="${bismark_bam/%.bam/.deduplication_report.txt}"


#########################


# exit if output exits already

if [ -s "$bismark_bam_dd_final" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 1
fi


#########################


# run deduplicate_bismark

module load bismark/0.16.3

echo " * BISMARK: $(readlink -f $(which deduplicate_bismark)) "
echo " * BAM IN: $bismark_bam "
echo " * BAM OUT ORIGINAL : $bismark_bam_dd_original "
echo " * REPORT ORIGINAL : $bismark_dd_report_original "
echo " * BAM OUT FINAL : $bismark_bam_dd_final "
echo " * REPORT FINAL : $bismark_dd_report_final "

if [ "$analysis_type" == "se" ] ; then
	bismark_flags="--single"
fi

if [ "$analysis_type" == "pe" ] ; then
	bismark_flags="--paired"
fi

bash_cmd="
deduplicate_bismark \
$bismark_flags \
--bam $bismark_bam
"
echo "CMD: $bash_cmd"
eval "$bash_cmd"

sleep 30


#########################


# check that output generated

if [ ! -s "$bismark_bam_dd_original" ] ; then
	echo -e "\n $script_name ERROR: REPORT $bismark_dd_report_original NOT GENERATED \n" >&2
	exit 1
fi

if [ ! -s "$bismark_dd_report_original" ] ; then
	echo -e "\n $script_name ERROR: REPORT $bismark_dd_report_original NOT GENERATED \n" >&2
	exit 1
fi


#########################


# put reports in report directory for bismark2report

bash_cmd="cp -v $bismark_dd_report_original $bismark_report_deduplication"
echo "CMD: $bash_cmd"
eval "$bash_cmd"


#########################


# clean up output name

bash_cmd="mv -v $bismark_bam_dd_original $bismark_bam_dd_final"
echo "CMD: $bash_cmd"
eval "$bash_cmd"

bash_cmd="mv -v $bismark_dd_report_original $bismark_dd_report_final"
echo "CMD: $bash_cmd"
eval "$bash_cmd"


#########################


# generate summary

reads_total=$(cat "$bismark_dd_report_final" | grep -m 1 "Total number of alignments analysed" | cut -f 2)
reads_duplicated=$(cat "$bismark_dd_report_final" | grep -m 1 "Total number duplicated alignments"  | cut -f 2 | cut -d ' ' -f 1)
reads_deduplicated=$(cat "$bismark_dd_report_final" | grep -m 1 "deduplicated leftover sequences" | cut -d ' ' -f 7)

reads_duplicated_pct=$(echo "(${reads_duplicated}/${reads_total})*100" | bc -l | cut -c 1-4)
reads_duplicated_pct="${reads_duplicated_pct}%"

# header
echo "#SAMPLE,ALIGNED PAIRS,DEDUPLICATED PAIRS,DUPLICATION RATE" > $summary_csv

# print the relevant numbers
echo "${sample},${reads_total},${reads_deduplicated},${reads_duplicated_pct}" >> "$summary_csv"

# paste -d ',' \
# <(echo "$sample") \
# <(cat "$bismark_dd_report_final" | grep -m 1 "Total number of alignments analysed" | cut -f 2) \
# <(cat "$bismark_dd_report_final" | grep -m 1 "deduplicated leftover sequences"     | cut -d ' ' -f 7) \
# <(cat "$bismark_dd_report_final" | grep -m 1 "Total number duplicated alignments"  | cut -f 2 | cut -d ' ' -f 1) \
# >> $summary_csv

sleep 30

# combine all sample summaries
cat ${summary_dir}/*.${segment_name}.csv | LC_ALL=C sort -t ',' -k1,1 | uniq \
> "${proj_dir}/summary.${segment_name}.csv"


#########################


# add sample and BAM to sample sheet
echo "${sample},${bismark_bam_dd_final}" >> "$samples_csv"

sleep 30

# sort and remove duplicates in place in sample sheet
LC_ALL=C sort -t ',' -k1,1 -u -o "$samples_csv" "$samples_csv"

sleep 30


#########################



# end
