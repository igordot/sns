#!/bin/bash


# run fastq_screen (more generic version of qc-fastqscreen.sh)


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

fastqscreen_conf="/gpfs/data/igorlab/software/fastq_screen_v0.13.0/fastq_screen.species.conf"

if [ ! -s "$fastqscreen_conf" ] ; then
	echo -e "\n $script_name ERROR: CONF $fastqscreen_conf DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# settings and files

fastqscreen_dir="${proj_dir}/species-fastqscreen"
mkdir -p "$fastqscreen_dir"

fastqscreen_txt="${fastqscreen_dir}/${fastq/*\//}_screen.txt"
fastqscreen_txt="${fastqscreen_txt/.fastq.gz/}"


#########################


# exit if output exists already

if [ -s "$fastqscreen_txt" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 1
fi


#########################


# fastq_screen
# ignore paired reads (in case of rna-seq, paired reads may be too far apart and will not align)

bowtie2_bin=$(cat "$fastqscreen_conf" | grep "^BOWTIE2" | head -1 | tr '[:space:]' '\t' | tr -s '\t' | cut -f 2)
fastqscreen_bin="/gpfs/data/igorlab/software/fastq_screen_v0.13.0/fastq_screen"

echo
echo " * fastq_screen: $fastqscreen_bin "
echo " * fastq_screen version: $($fastqscreen_bin --version) "
echo " * bowtie2: $bowtie2_bin "
echo " * bowtie2 version: $($bowtie2_bin --version 2>&1 | head -1) "
echo " * FASTQ: $fastq "
echo " * OUT TXT: $fastqscreen_txt "
echo

CMD="
$fastqscreen_bin \
--aligner bowtie2 \
--subset 1000000 \
--conf $fastqscreen_conf \
--outdir $fastqscreen_dir \
$fastq
"
echo "CMD: $CMD"
$CMD


#########################


# check that output generated

if [ ! -s "$fastqscreen_txt" ] ; then
	echo -e "\n $script_name ERROR: TXT $fastqscreen_txt NOT GENERATED \n" >&2
	exit 1
fi


#########################


# summary

# combine charts into a single pdf
combined_pdf="${proj_dir}/summary.fastqscreen.pdf"
rm -f "$combined_pdf"
convert "${fastqscreen_dir}/*screen.png" "$combined_pdf"


# combine charts into a single png

combined_png_3w="${proj_dir}/summary.fastqscreen.3w.png"
combined_png_4w="${proj_dir}/summary.fastqscreen.4w.png"

rm -f "$combined_png_3w"
rm -f "$combined_png_4w"

# -geometry +20+20 = 20px x and y padding
# -tile 3x = 3 images wide
montage -geometry +20+20 -tile 3x "${fastqscreen_dir}/*screen.png" "$combined_png_3w"
montage -geometry +20+20 -tile 4x "${fastqscreen_dir}/*screen.png" "$combined_png_4w"


#########################



# end
