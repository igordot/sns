#!/bin/bash


# run FastQ Screen


# script filename
script_path="${BASH_SOURCE[0]}"
script_name=$(basename "$script_path")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 4 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir sample_name num_threads FASTQ \n" >&2
	if [ $# -gt 0 ] ; then echo -e "\n ARGS: $* \n" >&2 ; fi
	exit 1
fi

# arguments
proj_dir=$1
sample=$2
threads=$3
fastq=$4


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

code_dir=$(dirname $(dirname "$script_path"))

fastqscreen_conf=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" REF-FASTQSCREEN);

if [ ! -s "$fastqscreen_conf" ] ; then
	echo -e "\n $script_name ERROR: CONF $fastqscreen_conf DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# settings and files

fastqscreen_dir="${proj_dir}/QC-fastqscreen"
mkdir -p "$fastqscreen_dir"

fastqscreen_txt="${fastqscreen_dir}/${fastq/*\//}_screen.txt"
fastqscreen_txt="${fastqscreen_txt/.fastq.gz/}"
fastqscreen_png="${fastqscreen_txt/_screen.txt/_screen.png}"

# unload all loaded modulefiles
module purge


#########################


# exit if output exists already

if [ -s "$fastqscreen_txt" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 0
fi


#########################


# fastq_screen
# ignore paired reads (in case of rna-seq, paired reads may be too far apart and will not align)

module add fastq_screen/0.13.0
# ImageMagick for "montage" for combining plots
module add imagemagick/7.0.8

bowtie2_bin=$(cat "$fastqscreen_conf" | grep "^BOWTIE2" | head -1 | tr '[:space:]' '\t' | tr -s '\t' | cut -f 2)

echo
echo " * fastq_screen: $(which fastq_screen) "
echo " * fastq_screen version: $(fastq_screen --version) "
echo " * bowtie2: $bowtie2_bin "
echo " * bowtie2 version: $($bowtie2_bin --version 2>&1 | head -1) "
echo " * fastq_screen conf: $fastqscreen_conf "
echo " * threads: $threads "
echo " * FASTQ: $fastq "
echo " * out TXT: $fastqscreen_txt "
echo " * out PNG: $fastqscreen_png "
echo

CMD="
fastq_screen \
--aligner bowtie2 \
--threads $threads \
--subset 1000000 \
--conf $fastqscreen_conf \
--outdir $fastqscreen_dir \
$fastq
"
echo "CMD: $CMD"
$CMD


#########################


# check that output generated

# check if TXT file is present
if [ ! -s "$fastqscreen_txt" ] ; then
	echo -e "\n $script_name ERROR: TXT $fastqscreen_txt NOT GENERATED \n" >&2
	exit 1
fi

# check if PNG file is present (if there was a problem with graphics devices)
if [ ! -s "$fastqscreen_png" ] ; then
	echo -e "\n $script_name ERROR: PNG $fastqscreen_png NOT GENERATED \n" >&2
	# delete TXT since something went wrong
	rm -fv "$fastqscreen_txt"
	exit 1
fi


#########################


# summary

# combine charts into a single pdf
combined_pdf="${proj_dir}/summary.fastqscreen.pdf"
rm -f "$combined_pdf"
convert "${fastqscreen_dir}/*screen.png" "$combined_pdf"


# combine charts into a single png

combined_png_2w="${proj_dir}/summary.fastqscreen.2w.png"
combined_png_4w="${proj_dir}/summary.fastqscreen.4w.png"

rm -f "$combined_png_2w"
rm -f "$combined_png_4w"

# -geometry +20+20 = 20px x and y padding
# -tile 3x = 3 images wide
montage -geometry +20+20 -tile 2x "${fastqscreen_dir}/*screen.png" "$combined_png_2w"
montage -geometry +20+20 -tile 4x "${fastqscreen_dir}/*screen.png" "$combined_png_4w"


#########################



# end
