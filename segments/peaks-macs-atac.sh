#!/bin/bash


# MACS peak calling for ATAC-seq


# script filename
script_name=$(basename "${BASH_SOURCE[0]}")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 4 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir sample_name BAM q_value \n" >&2
	exit 1
fi

# arguments
proj_dir=$1
sample=$2
bam=$3
qvalue=$4


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

code_dir=$(dirname "$(dirname "${BASH_SOURCE[0]}")")

genome_dir=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" GENOME-DIR);

if [ ! -d "$genome_dir" ] ; then
	echo -e "\n $script_name ERROR: GENOME DIR $genome_dir DOES NOT EXIST \n" >&2
	exit 1
fi

chrom_sizes=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" REF-CHROMSIZES);

if [ ! -s "$chrom_sizes" ] ; then
	echo -e "\n $script_name ERROR: CHROM SIZES $chrom_sizes DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# MACS-stlye genome (mm,hs,dm,ce)

# keep first two characters of build name
genome_build=$(basename "$genome_dir")
macs_genome="${genome_build:0:2}"

# fix if hgXX
if [ "$macs_genome" == "hg" ] ; then
	macs_genome="hs"
fi


#########################


# settings and files

segment_name="${segment_name}-q-${qvalue}"

summary_dir="${proj_dir}/summary"
mkdir -p "$summary_dir"
summary_csv="${summary_dir}/${sample}.${segment_name}.csv"

macs_dir="${proj_dir}/peaks-macs-atac-q-${qvalue}"
mkdir -p "$macs_dir"
peaks_xls="${macs_dir}/${sample}_peaks.xls"
peaks_narrow="${macs_dir}/${sample}_peaks.narrowPeak"
peaks_bed="${macs_dir}/${sample}.bed"
macs_bdg_treat="${macs_dir}/${sample}_treat_pileup.bdg"
macs_bdg_treat_sort="${macs_dir}/${sample}_treat_pileup.sort.bdg"
macs_bdg_control="${macs_dir}/${sample}_control_lambda.bdg"

bigwig_dir="${proj_dir}/BIGWIG"
mkdir -p "$bigwig_dir"
macs_bw="${bigwig_dir}/${sample}.macs.bw"


#########################


# exit if output exits already

if [ -s "$peaks_xls" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 1
fi


#########################


# MACS

# MACS is part of python/2.7.3 module
module unload python
module load python/2.7.3

echo " * MACS: $(readlink -f $(which macs2)) "
echo " * MACS version: $(macs2 --version 2>&1) "
echo " * BAM: $bam "
echo " * MACS genome: $macs_genome "
echo " * MACS dir: $macs_dir "

cd "$macs_dir" || exit 1

bash_cmd="
macs2 callpeak \
--verbose 2 \
--nomodel --shift -100 --extsize 200 \
--bdg --SPMR \
--keep-dup all \
--qvalue $qvalue \
--gsize $macs_genome \
--name $sample \
--treatment $bam \
--outdir $macs_dir
"
echo -e "\n CMD: $bash_cmd \n"
$bash_cmd


#########################


# check that output generated

if [ ! -s "$peaks_xls" ] ; then
	echo -e "\n $script_name ERROR: XLS $peaks_xls NOT GENERATED \n" >&2
	exit 1
fi

if [ ! -s "$peaks_narrow" ] ; then
	echo -e "\n $script_name ERROR: PEAKS $peaks_narrow NOT GENERATED \n" >&2
	exit 1
fi


#########################


# generate bed file

bash_cmd="cut -f 1,2,3,4,7 $peaks_narrow | LC_ALL=C sort -k1,1 -k2,2n > $peaks_bed"
echo -e "\n CMD: $bash_cmd \n"
eval "$bash_cmd"


#########################


# generate bigwig file if not already generated

if [ ! -s "$macs_bdg_treat" ] ; then
	echo -e "\n $script_name ERROR: BEDGRAPH $macs_bdg_treat NOT GENERATED \n" >&2
	exit 1
fi

module load kentutils/329

if [ ! -s "$macs_bw" ] ; then

	echo " * bedGraphToBigWig: $(readlink -f $(which bedGraphToBigWig)) "
	echo " * MACS bedGraph original: $macs_bdg_treat "
	echo " * MACS bedGraph sorted: $macs_bdg_treat_sort "
	echo " * MACS bigWig: $macs_bw "

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

peaks_num=$(cat "$peaks_narrow" | wc -l)
echo "total peaks: $peaks_num"

# header for summary file
echo "#SAMPLE,PEAKS q ${qvalue}" > "$summary_csv"

# summarize log file
echo "${sample},${peaks_num}" >> "$summary_csv"

sleep 30

# combine all sample summaries
cat ${summary_dir}/*.${segment_name}.csv | LC_ALL=C sort -t ',' -k1,1 | uniq > "${proj_dir}/summary.${segment_name}.csv"


#########################



# end
