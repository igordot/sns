#!/bin/bash


# run Picard CollectRnaSeqMetrics


# script filename
script_name=$(basename "${BASH_SOURCE[0]}")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 3 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name [project dir] [sample] [BAM] \n" >&2
	exit 1
fi

# arguments
proj_dir=$1
sample=$2
bam=$3

# strand:
# fwd | transcript             | cufflinks "fr-secondstrand" | htseq "yes"     | picard "FIRST_READ"
# rev | rev comp of transcript | cufflinks "fr-firststrand"  | htseq "reverse" | picard "SECOND_READ"


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

CODE_DIR=$(dirname "$(dirname "${BASH_SOURCE[0]}")")
REFFLAT=$(bash ${CODE_DIR}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" REF-REFFLAT)
RRNA_INTERVAL_LIST=$(bash ${CODE_DIR}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" REF-RRNAINTERVALLIST)

if [ ! -s "$REFFLAT" ] || [ ! "$REFFLAT" ] ; then
	echo -e "\n $script_name ERROR: REFFLAT $REFFLAT DOES NOT EXIST \n" >&2
	exit 1
fi

if [ ! -s "$RRNA_INTERVAL_LIST" ] || [ ! "$RRNA_INTERVAL_LIST" ] ; then
	echo -e "\n $script_name ERROR: RRNA INTERVAL LIST $RRNA_INTERVAL_LIST DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# output

metrics_dir="${proj_dir}/QC-RnaSeqMetrics"
mkdir -p "$metrics_dir"
metrics_txt="${metrics_dir}/${sample}.txt"
metrics_pdf="${metrics_dir}/${sample}.pdf"

summary_dir="${proj_dir}/summary"
summary_csv="${summary_dir}/${sample}.${segment_name}.csv"

# # strand flag
# if [ "$STRAND" == "unstr" ] ; then
# 	STRAND_FLAG="0"
# elif [ "$STRAND" == "fwd" ] ; then
# 	STRAND_FLAG="1"
# elif [ "$STRAND" == "rev" ] ; then
# 	STRAND_FLAG="2"
# else
# 	echo -e "\n $script_name ERROR: incorrect strand selected \n" >&2
# 	exit 1
# fi


#########################


# exit if output exits already

if [ -s "$metrics_txt" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 1
fi


#########################


module unload java
module load java/1.7
module load picard-tools/1.88




# STRAND_SPECIFICITY {NONE, FIRST_READ_TRANSCRIPTION_STRAND, SECOND_READ_TRANSCRIPTION_STRAND}
# STRAND_SPECIFICITY="NONE"


echo " * CollectRnaSeqMetrics: ${PICARD_ROOT}/CollectRnaSeqMetrics.jar "
echo " * BAM: $bam "
echo " * REF_FLAT: $REFFLAT "
echo " * RIBOSOMAL_INTERVALS: $RRNA_INTERVAL_LIST "
echo " * OUT TXT: $metrics_txt "
echo " * OUT PDF: $metrics_pdf "

# run picard

PICARD_CMD="java -Xms8G -Xmx8G -jar ${PICARD_ROOT}/CollectRnaSeqMetrics.jar \
VERBOSITY=WARNING QUIET=true VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=2500000 \
REF_FLAT=${REFFLAT} \
INPUT=${bam}"

CMD="$PICARD_CMD \
STRAND_SPECIFICITY=NONE \
RIBOSOMAL_INTERVALS=${RRNA_INTERVAL_LIST} \
CHART_OUTPUT=${metrics_pdf} \
OUTPUT=${metrics_txt}"
echo "CMD: $CMD"
$CMD

CMD="$PICARD_CMD \
STRAND_SPECIFICITY=FIRST_READ_TRANSCRIPTION_STRAND \
OUTPUT=${metrics_txt}.1READ"
echo "CMD: $CMD"
$CMD

CMD="$PICARD_CMD \
STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
OUTPUT=${metrics_txt}.2READ"
echo "CMD: $CMD"
$CMD


#########################


# check that output generated

if [ ! -s "$metrics_txt" ] ; then
	echo -e "\n $script_name ERROR: METRICS $metrics_txt NOT GENERATED \n" >&2
	exit 1
fi


#########################


# generate a summary file

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
# 11 PCT_RIBOSOMAL_BASES
# 12 PCT_CODING_BASES
# 13 PCT_UTR_BASES
# 14 PCT_INTRONIC_BASES
# 15 PCT_INTERGENIC_BASES
# 16 PCT_MRNA_BASES
# 17 PCT_USABLE_BASES
# 18 PCT_CORRECT_STRAND_READS
# 19 MEDIAN_CV_COVERAGE
# 20 MEDIAN_5PRIME_BIAS
# 21 MEDIAN_3PRIME_BIAS
# 22 MEDIAN_5PRIME_TO_3PRIME_BIAS

paste \
<(echo -e "#SAMPLE\n${sample}") \
<(cat "${metrics_txt}" | grep -A 1 "PF_ALIGNED_BASES" | cut -f 2,11,12,13,14,15,22) \
<(cat "${metrics_txt}.1READ" | grep -A 1 "PF_ALIGNED_BASES" | cut -f 18 | sed 's/PCT_CORRECT_STRAND_READS/FWD_STRAND/') \
<(cat "${metrics_txt}.2READ" | grep -A 1 "PF_ALIGNED_BASES" | cut -f 18 | sed 's/PCT_CORRECT_STRAND_READS/REV_STRAND/') \
| tr '\t' ',' \
> $summary_csv

rm -f "${metrics_txt}.1READ"
rm -f "${metrics_txt}.2READ"


# combine charts into a single pdf

combined_pdf=${proj_dir}/summary.${segment_name}.pdf

rm -f "$combined_pdf"

# do not use quotes in filenames for ghostscript
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=${combined_pdf} ${metrics_dir}/*.pdf

# combine charts into a single png

combined_png_4w=${proj_dir}/summary.${segment_name}.4w.png
combined_png_5w=${proj_dir}/summary.${segment_name}.5w.png

rm -f "$combined_png_4w"
rm -f "$combined_png_5w"

# -geometry +20+20 = 20px x and y padding
# -tile 4x = 4 images wide
montage -geometry +20+20 -tile 4x "${metrics_dir}/*.pdf" "$combined_png_4w"
montage -geometry +20+20 -tile 5x "${metrics_dir}/*.pdf" "$combined_png_5w"


# combine all sample summaries
cat ${summary_dir}/*.${segment_name}.csv| LC_ALL=C sort -t ',' -k1,1 | uniq > "${proj_dir}/summary.${segment_name}.csv"


#########################



# end
