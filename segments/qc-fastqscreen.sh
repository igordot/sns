#!/bin/bash


# run fastq_screen


# script filename
script_name=$(basename "${BASH_SOURCE[0]}")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 3 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name [project dir] [sample] [FASTQ] \n" >&2
	exit 1
fi

# arguments
proj_dir=$1
sample=$2
fastq=$3


#########################


# check that inputs exist

if [ ! -d "$proj_dir" ] ; then
	echo -e "\n $script_name ERROR: DIR $proj_dir DOES NOT EXIST \n" >&2
	exit 1
fi

if [ ! -s "$fastq" ] ; then
	echo -e "\n $script_name ERROR: FASTQ $fastq DOES NOT EXIST \n" >&2
	exit 1
fi

CODE_DIR=$(dirname "$(dirname "${BASH_SOURCE[0]}")")
FASTQSCREEN_CONF=$(bash ${CODE_DIR}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" REF-FASTQSCREEN);

if [ ! -s "$FASTQSCREEN_CONF" ] ; then
	echo -e "\n $script_name ERROR: CONF $FASTQSCREEN_CONF DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# settings and files

FASTQ_SCREEN_DIR="${proj_dir}/QC-fastqscreen"
mkdir -p "$FASTQ_SCREEN_DIR"

FASTQ_SCREEN_TXT="${FASTQ_SCREEN_DIR}/${fastq/*\//}_screen.txt"
FASTQ_SCREEN_TXT="${FASTQ_SCREEN_TXT/.fastq.gz/}"


#########################


# exit if output exits already

if [ -s "$FASTQ_SCREEN_TXT" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 1
fi


#########################


# fastq_screen
# ignore paired reads (in case of rna-seq, paired reads may be too far apart and will not align)

bowtie2_bin=$(cat "$FASTQSCREEN_CONF" | grep "^BOWTIE2" | head -1 | tr '[:space:]' '\t' | tr -s '\t' | cut -f 2)
fastqscreen_bin="/ifs/home/id460/software/fastq_screen_v0.5.2/fastq_screen"

echo " * fastq_screen: $fastqscreen_bin "
echo " * fastq_screen version: $($fastqscreen_bin --version) "
echo " * bowtie2: $bowtie2_bin "
echo " * bowtie2 version: $($bowtie2_bin --version 2>&1 | head -1) "
echo " * FASTQ: $fastq "
echo " * OUT TXT: $FASTQ_SCREEN_TXT "

CMD="
$fastqscreen_bin \
--aligner bowtie2 \
--subset 1000000 \
--conf $FASTQSCREEN_CONF \
--outdir $FASTQ_SCREEN_DIR \
$fastq
"
echo "CMD: $CMD"
$CMD


#########################


# check that output generated

if [ ! -s "$FASTQ_SCREEN_TXT" ] ; then
	echo -e "\n $script_name ERROR: TXT $FASTQ_SCREEN_TXT NOT GENERATED \n" >&2
	exit 1
fi


#########################


# summary

# combine charts into a single pdf
COMBINED_PDF="${proj_dir}/summary.fastqscreen.pdf"
rm -f "$COMBINED_PDF"
convert "${FASTQ_SCREEN_DIR}/*screen.png" "$COMBINED_PDF"


# combine charts into a single png

COMBINED_PNG_3W="${proj_dir}/summary.fastqscreen.3w.png"
COMBINED_PNG_4W="${proj_dir}/summary.fastqscreen.4w.png"

rm -f "$COMBINED_PNG_3W"
rm -f "$COMBINED_PNG_4W"

# -geometry +20+20 = 20px x and y padding
# -tile 3x = 3 images wide
montage -geometry +20+20 -tile 3x "${FASTQ_SCREEN_DIR}/*screen.png" "$COMBINED_PNG_3W"
montage -geometry +20+20 -tile 4x "${FASTQ_SCREEN_DIR}/*screen.png" "$COMBINED_PNG_4W"


#########################



# end
