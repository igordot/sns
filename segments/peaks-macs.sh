#!/bin/bash


# MACS peak calling


# script filename
script_path="${BASH_SOURCE[0]}"
script_name=$(basename "$script_path")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ $# -lt 5 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir peak_type peak_q sample_name treatment_bam [control_bam] \n" >&2
	exit 1
fi

# arguments
proj_dir=$(readlink -f "$1")
peak_type=$2
peak_q=$3
sample=$4
bam_treat=$5
bam_control=$6


#########################


# check that inputs exist

if [ ! -d "$proj_dir" ] ; then
	echo -e "\n $script_name ERROR: PROJ DIR $proj_dir DOES NOT EXIST \n" >&2
	exit 1
fi

if [ ! -s "$bam_treat" ] ; then
	echo -e "\n $script_name ERROR: BAM $bam_treat DOES NOT EXIST \n" >&2
	exit 1
fi

code_dir=$(dirname $(dirname "$script_path"))

genome_dir=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" GENOME-DIR);

if [ ! -d "$genome_dir" ] ; then
	echo -e "\n $script_name ERROR: GENOME DIR $genome_dir DOES NOT EXIST \n" >&2
	exit 1
fi

chrom_sizes=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" REF-CHROMSIZES);

if [ ! -s "$chrom_sizes" ] ; then
	echo -e "\n $script_name ERROR: CHROM SIZES FILE $chrom_sizes DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# settings and files

segment_name="${segment_name}-${peak_type}-q-${peak_q}"
if [ -z "$bam_control" ] ; then
	segment_name="${segment_name}-nocontrol"
fi

summary_dir="${proj_dir}/summary"
mkdir -p "$summary_dir"
summary_csv="${summary_dir}/${sample}.${segment_name}.csv"

macs_dir="${proj_dir}/${segment_name}"
mkdir -p "$macs_dir"
peaks_xls="${macs_dir}/${sample}_peaks.xls"
peaks_file="${macs_dir}/${sample}_peaks.narrowPeak"
peaks_bed="${macs_dir}/${sample}.bed"
model_r="${macs_dir}/${sample}_model.r"
macs_bdg_treat="${macs_dir}/${sample}_treat_pileup.bdg"
macs_bdg_treat_sort="${macs_dir}/${sample}_treat_pileup.sort.bdg"
macs_bdg_control="${macs_dir}/${sample}_control_lambda.bdg"

bigwig_dir="${proj_dir}/BIGWIG"
mkdir -p "$bigwig_dir"
macs_bw="${bigwig_dir}/${sample}.macs.bw"

# unload all loaded modulefiles
module purge
module add default-environment


#########################


# exit if output exits already

if [ -s "$peaks_xls" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 1
fi


#########################


# MACS parameters

# keep first two characters of build name
genome_build=$(basename "$genome_dir")
macs_genome="${genome_build:0:2}"

# fix if hgXX
if [ "$macs_genome" == "hg" ] ; then
	macs_genome="hs"
fi

# adjust for control sample
if [ -n "$bam_control" ] ; then
	bam_param="--treatment $bam_treat --control $bam_control"
else
	bam_param="--treatment $bam_treat"
fi

# adjust parameters based on run type
if [ "$peak_type" == "atac" ] ; then
	analysis_params="--nomodel --shift -100 --extsize 200 --qvalue $peak_q"
elif [ "$peak_type" == "narrow" ] ; then
	analysis_params="--qvalue $peak_q"
elif [ "$peak_type" == "broad" ] ; then
	analysis_params="--broad --broad-cutoff $peak_q"
	peaks_file="${macs_dir}/${sample}_peaks.broadPeak"
else
	echo -e "\n $script_name ERROR: unknown peak type $peak_type \n" >&2
	exit 1
fi


#########################


# MACS

# MACS is part of python/cpu/2.7.15 module
module add python/cpu/2.7.15

echo
echo " * MACS: $(readlink -f $(which macs2)) "
echo " * MACS version: $(macs2 --version 2>&1) "
echo " * BAM treatment: $bam_treat "
echo " * BAM control: $bam_control "
echo " * MACS genome: $macs_genome "
echo " * MACS dir: $macs_dir "
echo " * peak type: $peak_type "
echo " * peak cutoff: $peak_q "
echo

cd "$macs_dir" || exit 1

bash_cmd="
macs2 callpeak \
--verbose 2 \
--bdg --SPMR \
--keep-dup all \
$analysis_params \
--gsize $macs_genome \
--name $sample \
$bam_param \
--outdir $macs_dir
"
echo -e "\n CMD: $bash_cmd \n"
$bash_cmd


#########################


# check that output generated

if [ ! -s "$peaks_file" ] ; then
	echo -e "\n $script_name ERROR: PEAKS FILE $peaks_file NOT GENERATED \n" >&2
	exit 1
fi

if [ ! -s "$peaks_xls" ] ; then
	echo -e "\n $script_name ERROR: XLS $peaks_xls NOT GENERATED \n" >&2
	exit 1
fi

if [ ! -s "$model_r" ] ; then
	echo -e "\n $script_name ERROR: MODEL R SCRIPT $model_r NOT GENERATED \n" >&2
	exit 1
fi


#########################


# generate an image about the model based on the data

module add r/3.6.1

echo
echo " * R: $(readlink -f $(which R)) "
echo " * R version: $(R --version | head -1) "
echo " * Rscript: $(readlink -f $(which Rscript)) "
echo " * Rscript version: $(Rscript --version 2>&1) "
echo

bash_cmd="Rscript $model_r"
echo -e "\n CMD: $bash_cmd \n"
$bash_cmd

sleep 5


#########################


# generate bed file

bash_cmd="cut -f 1,2,3,4,7 $peaks_file | LC_ALL=C sort -k1,1 -k2,2n > $peaks_bed"
echo -e "\n CMD: $bash_cmd \n"
eval "$bash_cmd"

sleep 5


#########################


# generate bigwig file if not already generated

if [ ! -s "$macs_bdg_treat" ] ; then
	echo -e "\n $script_name ERROR: BEDGRAPH $macs_bdg_treat NOT GENERATED \n" >&2
	exit 1
fi

module add ucscutils/374

if [ ! -s "$macs_bw" ] ; then

	echo
	echo " * bedGraphToBigWig: $(readlink -f $(which bedGraphToBigWig)) "
	echo " * MACS bedGraph original: $macs_bdg_treat "
	echo " * MACS bedGraph sorted: $macs_bdg_treat_sort "
	echo " * MACS bigWig: $macs_bw "
	echo

	bdg_sort_cmd="cat $macs_bdg_treat | LC_ALL=C sort -k1,1 -k2,2n > $macs_bdg_treat_sort"
	echo -e "\n CMD: $bdg_sort_cmd \n"
	eval "$bdg_sort_cmd"

	bw_cmd="bedGraphToBigWig $macs_bdg_treat_sort $chrom_sizes $macs_bw"
	echo -e "\n CMD: $bw_cmd \n"
	eval "$bw_cmd"

	rm -fv "$macs_bdg_treat_sort"

fi

rm -fv "$macs_bdg_treat"
rm -fv "$macs_bdg_control"


#########################


# generate summary

peaks_num=$(cat "$peaks_bed" | wc -l)
echo "total peaks: $peaks_num"

# header for summary file
echo "#SAMPLE,PEAKS q ${peak_q}" > "$summary_csv"

# summarize log file
echo "${sample},${peaks_num}" >> "$summary_csv"

sleep 5

# combine all sample summaries
cat ${summary_dir}/*.${segment_name}.csv | LC_ALL=C sort -t ',' -k1,1 | uniq > "${proj_dir}/summary.${segment_name}.csv"


#########################



# end
