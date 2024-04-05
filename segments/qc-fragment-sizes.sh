#!/bin/bash


# get fragment size distribution


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
proj_dir=$(readlink -f "$1")
sample=$2
bam=$3


#########################


# settings and files

summary_dir="${proj_dir}/summary"
mkdir -p "$summary_dir"
summary_csv="${summary_dir}/${sample}.${segment_name}.csv"

frag_sizes_dir="${proj_dir}/QC-fragment-sizes"
mkdir -p "$frag_sizes_dir"
frag_sizes_png="${frag_sizes_dir}/${sample}.sizes.600.png"
frag_sizes_csv="${frag_sizes_dir}/${sample}.sizes.csv"
frag_sizes_stats_csv="${frag_sizes_dir}/${sample}.stats.csv"

# unload all loaded modulefiles
module purge


#########################


# exit if output exits already

if [ -s "$frag_sizes_csv" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 0
fi


#########################


# check that inputs exist

if [ ! -d "$proj_dir" ] ; then
	echo -e "\n $script_name ERROR: PROJ DIR $proj_dir DOES NOT EXIST \n" >&2
	exit 1
fi

if [ ! -s "$bam" ] ; then
	echo -e "\n $script_name ERROR: BAM $bam DOES NOT EXIST \n" >&2
	exit 1
fi

code_dir=$(dirname $(dirname "$script_path"))


#########################


# calculate and plot fragment size distribution

module add r/4.1.2

echo
echo " * R: $(readlink -f $(which R)) "
echo " * R version: $(R --version | head -1) "
echo " * Rscript: $(readlink -f $(which Rscript)) "
echo " * Rscript version: $(Rscript --version 2>&1) "
echo

# test R
Rscript --vanilla "${code_dir}/scripts/test-package.R" getopt
Rscript --vanilla "${code_dir}/scripts/test-package.R" optparse

# navigate to frag sizes dir
cd "$frag_sizes_dir" || exit 1

# run fragment sizes script
bash_cmd="Rscript --vanilla ${code_dir}/scripts/fragment-sizes.R $sample $bam"
echo "CMD: $bash_cmd"
($bash_cmd)


#########################


# check that output generated

# check if CSV file is present
if [ ! -s "$frag_sizes_csv" ] ; then
	echo -e "\n $script_name ERROR: FILE $frag_sizes_csv NOT GENERATED \n" >&2
	exit 1
fi

# check if PNG file is present (if there was a problem with graphics devices)
if [ ! -s "$frag_sizes_png" ] ; then
	echo -e "\n $script_name ERROR: FILE $frag_sizes_png NOT GENERATED \n" >&2
	# delete CSV since something went wrong
	rm -fv "$frag_sizes_csv"
	exit 1
fi


#########################


# summary

# "montage" is part of ImageMagick
module add imagemagick/7.0.8

# combine charts into a single png

combined_png_2w=${proj_dir}/summary.${segment_name}.2w.png
combined_png_4w=${proj_dir}/summary.${segment_name}.4w.png

rm -f "$combined_png_2w"
rm -f "$combined_png_4w"

# -geometry +20+20 = 20px x and y padding
# -tile 4x = 4 images wide
montage -geometry +20+20 -tile 2x "${frag_sizes_dir}/*.600.png" "$combined_png_2w"
montage -geometry +20+20 -tile 4x "${frag_sizes_dir}/*.600.png" "$combined_png_4w"

# header for summary file
echo "#SAMPLE,MEAN FRAGMENT,MEDIAN FRAGMENT,SD FRAGMENT" > "$summary_csv"

# summarize log file
grep -v "^SAMPLE" "$frag_sizes_stats_csv" >> "$summary_csv"

sleep 5

# combine all sample summaries
cat ${summary_dir}/*.${segment_name}.csv | LC_ALL=C sort -t ',' -k1,1 | uniq > "${proj_dir}/summary.${segment_name}.csv"


#########################



# end
