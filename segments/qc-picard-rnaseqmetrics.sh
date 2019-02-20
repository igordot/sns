#!/bin/bash


# run Picard CollectRnaSeqMetrics


# script filename
script_path="${BASH_SOURCE[0]}"
script_name=$(basename "$script_path")
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
proj_dir=$1
sample=$2
bam=$3


#########################


# settings and files

metrics_dir="${proj_dir}/QC-RnaSeqMetrics"
mkdir -p "$metrics_dir"
metrics_txt="${metrics_dir}/${sample}.txt"
metrics_pdf="${metrics_dir}/${sample}.pdf"

summary_dir="${proj_dir}/summary"
summary_csv="${summary_dir}/${sample}.${segment_name}.csv"

# unload all loaded modulefiles
module purge
module add default-environment


#########################


# exit if output exits already

if [ -s "$metrics_txt" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 1
fi


#########################


# check that inputs exist

if [ ! -d "$proj_dir" ] || [ ! "$proj_dir" ] ; then
	echo -e "\n $script_name ERROR: DIR $proj_dir DOES NOT EXIST \n" >&2
	exit 1
fi

if [ ! -s "$bam" ] || [ ! "$bam" ] ; then
	echo -e "\n $script_name ERROR: BAM $bam DOES NOT EXIST \n" >&2
	exit 1
fi

code_dir=$(dirname $(dirname "$script_path"))

refflat=$(bash "${code_dir}/scripts/get-set-setting.sh" "${proj_dir}/settings.txt" REF-REFFLAT)

rrna_interval_list=$(bash "${code_dir}/scripts/get-set-setting.sh" "${proj_dir}/settings.txt" REF-RRNAINTERVALLIST)

if [ ! -s "$refflat" ] || [ ! "$refflat" ] ; then
	echo -e "\n $script_name ERROR: REFFLAT $refflat DOES NOT EXIST \n" >&2
	exit 1
fi

if [ ! -s "$rrna_interval_list" ] || [ ! "$rrna_interval_list" ] ; then
	echo -e "\n $script_name ERROR: RRNA INTERVAL LIST $rrna_interval_list DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# run Picard CollectRnaSeqMetrics

module add picard-tools/2.18.20
# Picard requires R for some plotting
module add r/3.5.1
# ImageMagick for "montage" for combining plots
module add imagemagick/7.0.8

picard_jar="${PICARD_ROOT}/libs/picard.jar"

# check if the picard jar file is present
if [ ! -s "$picard_jar" ] ; then
	echo -e "\n $script_name ERROR: FILE $picard_jar DOES NOT EXIST \n" >&2
	exit 1
fi

# strand                       | cufflinks       | htseq   | picard      |
# fwd (transcript)             | fr-secondstrand | yes     | FIRST_READ  |
# rev (rev comp of transcript) | fr-firststrand  | reverse | SECOND_READ |

echo
echo " * Picard: $picard_jar "
echo " * BAM: $bam "
echo " * refFlat: $refflat "
echo " * ribosomal intervals: $rrna_interval_list "
echo " * out TXT: $metrics_txt "
echo " * out PDF: $metrics_pdf "
echo

picard_base_cmd="java -Xms16G -Xmx16G -jar $picard_jar CollectRnaSeqMetrics \
VERBOSITY=WARNING QUIET=true VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=2500000 \
REF_FLAT=${refflat} \
INPUT=${bam}"

picard_unstranded_cmd="$picard_base_cmd \
STRAND_SPECIFICITY=NONE \
RIBOSOMAL_INTERVALS=${rrna_interval_list} \
CHART_OUTPUT=${metrics_pdf} \
OUTPUT=${metrics_txt}"
echo "CMD: $picard_unstranded_cmd"
$picard_unstranded_cmd


#########################


# check that output generated

if [ ! -s "$metrics_txt" ] ; then
	echo -e "\n $script_name ERROR: METRICS $metrics_txt NOT GENERATED \n" >&2
	exit 1
fi


#########################


# generate a summary file

# RnaSeqMetrics output
# https://broadinstitute.github.io/picard/picard-metric-definitions.html#RnaSeqMetrics

# 01 PF_BASES
# 02 PF_ALIGNED_BASES
# 03 RIBOSOMAL_BASES
# 04 CODING_BASES
# 05 UTR_BASES
# 06 INTRONIC_BASES
# 07 INTERGENIC_BASES
# 08 IGNORED_READS
# 09 CORRECT_STRAND_READS
# 10 INCORRECT_STRAND_READS
# 11 NUM_R1_TRANSCRIPT_STRAND_READS
# 12 NUM_R2_TRANSCRIPT_STRAND_READS
# 13 NUM_UNEXPLAINED_READS
# 14 PCT_R1_TRANSCRIPT_STRAND_READS
# 15 PCT_R2_TRANSCRIPT_STRAND_READS
# 16 PCT_RIBOSOMAL_BASES
# 17 PCT_CODING_BASES
# 18 PCT_UTR_BASES
# 19 PCT_INTRONIC_BASES
# 20 PCT_INTERGENIC_BASES
# 21 PCT_MRNA_BASES
# 22 PCT_USABLE_BASES
# 23 PCT_CORRECT_STRAND_READS
# 24 MEDIAN_CV_COVERAGE
# 25 MEDIAN_5PRIME_BIAS
# 26 MEDIAN_3PRIME_BIAS
# 27 MEDIAN_5PRIME_TO_3PRIME_BIAS
# 28 SAMPLE
# 29 LIBRARY
# 30 READ_GROUP

# subset for specific metrics
# splitting by read may not be needed with Picard 2.18
paste \
<(echo -e "#SAMPLE\n${sample}") \
<(cat "${metrics_txt}" | grep -A 1 "PF_ALIGNED_BASES" | cut -f 2,14,15,17,18,19,20,21,27) \
| tr '\t' ',' \
> $summary_csv


# combine charts into a single pdf

combined_pdf="${proj_dir}/summary.${segment_name}.pdf"

rm -f "$combined_pdf"

# do not use quotes in filenames for ghostscript
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=${combined_pdf} ${metrics_dir}/*.pdf

# combine charts into a single png

combined_png_3w="${proj_dir}/summary.${segment_name}.3w.png"
combined_png_5w="${proj_dir}/summary.${segment_name}.5w.png"

rm -f "$combined_png_3w"
rm -f "$combined_png_5w"

# -geometry +20+20 = 20px x and y padding
# -tile 4x = 4 images wide
montage -geometry +20+20 -tile 3x "${metrics_dir}/*.pdf" "$combined_png_3w"
montage -geometry +20+20 -tile 5x "${metrics_dir}/*.pdf" "$combined_png_5w"


# combine all sample summaries
cat ${summary_dir}/*.${segment_name}.csv| LC_ALL=C sort -t ',' -k1,1 | uniq > "${proj_dir}/summary.${segment_name}.csv"


#########################



# end
