#!/bin/bash


# HMMRATAC peak calling for ATAC-seq data


# script filename
script_path="${BASH_SOURCE[0]}"
script_name=$(basename "$script_path")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 3 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir sample_name bam \n" >&2
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

macs_peaks_dir="${proj_dir}/peaks-MACS3-HMMRATAC"
mkdir -p "$macs_peaks_dir"
peaks_bed="${macs_peaks_dir}/${sample}.peaks.bed"

macs_logs_dir="${proj_dir}/logs-${segment_name}"
mkdir -p "$macs_logs_dir"
macs_log_txt="${macs_logs_dir}/${sample}.macs3.txt"
peaks_cutoff_tsv="${macs_logs_dir}/${sample}_cutoff_analysis.tsv"
peaks_model_json="${macs_logs_dir}/${sample}_model.json"
# gappedPeak in MACS 3.0.1, narrowPeak in MACS 3.0.2
peaks_file="${macs_logs_dir}/${sample}_accessible_regions.narrowPeak"

# unload all loaded modulefiles
module purge


#########################


# exit if output exits already

if [ -s "$peaks_bed" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 1
fi


#########################


# check that inputs exist

if [ ! -d "$proj_dir" ] ; then
	echo -e "\n $script_name ERROR: proj dir $proj_dir does not exist \n" >&2
	exit 1
fi

if [ ! -s "$bam" ] ; then
	echo -e "\n $script_name ERROR: BAM $bam does not exist \n" >&2
	exit 1
fi

code_dir=$(dirname $(dirname "$script_path"))

genome_dir=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" GENOME-DIR);

if [ ! -d "$genome_dir" ] ; then
	echo -e "\n $script_name ERROR: genome dir $genome_dir does not exist \n" >&2
	exit 1
fi

blacklist=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" REF-BLACKLIST);

if [ ! -s "$blacklist" ] ; then
	echo -e "\n $script_name ERROR: blacklist $blacklist does not exist \n" >&2
	exit 1
fi


#########################


# MACS

# MACS3 is part of condaenvs/2023/macs3 module
module add condaenvs/2023/macs3

echo
echo " * MACS path: $(readlink -f $(which macs3)) "
echo " * MACS version: $(macs3 --version 2>&1) "
echo " * Python path: $(readlink -f $(which python)) "
echo " * Python version: $(python --version 2>&1) "
echo " * blacklist: $blacklist "
echo " * BAM: $bam "
echo " * MACS logs dir: $macs_logs_dir "
echo " * MACS peaks dir: $macs_peaks_dir "
echo

cd "$macs_logs_dir" || exit 1

bash_cmd="
macs3 hmmratac \
--input $bam \
--name $sample \
--blacklist $blacklist \
--outdir $macs_logs_dir \
--hmm-type poisson \
2> $macs_log_txt \
"
echo -e "\n CMD: $bash_cmd \n"
eval "$bash_cmd"


#########################


# check that output generated

if [ ! -s "$peaks_cutoff_tsv" ] ; then
	echo -e "\n $script_name ERROR: cutoff analysis $peaks_cutoff_tsv not generated \n" >&2
	exit 1
fi

if [ ! -s "$peaks_model_json" ] ; then
	echo -e "\n $script_name ERROR: model JSON $peaks_model_json not generated \n" >&2
	exit 1
fi

if [ ! -s "$peaks_file" ] ; then
	echo -e "\n $script_name ERROR: narrowPeak $peaks_file not generated \n" >&2
	exit 1
fi


#########################


# generate a clean BED file from a narrowPeak file

# narrowPeak format:
# chrom start end name score strand signalValue pValue qValue peak

# hmmratac narrowPeak format (accessible regions or HMM open state):
# 1: chromosome name
# 2: start position of the accessible region
# 3: end position of the accessible region
# 4: peak name
# 5: peak score - maximum fold change (signal/average signal) within the peak
# 6: not used
# 7: not used
# 8: not used
# 9: peak summit position - relative position from the start to the peak summit (maximum fold change)

bash_cmd="
cut -f 1,2,3,4 $peaks_file \
| grep -v 'type=narrowPeak' \
| cut -f 1-6 \
| LC_ALL=C sort -k1,1 -k2,2n \
> $peaks_bed
"
echo -e "\n CMD: $bash_cmd \n"
eval "$bash_cmd"


#########################


# check that output generated

if [ ! -s "$peaks_bed" ] ; then
	echo -e "\n $script_name ERROR: peaks $peaks_bed not generated \n" >&2
	exit 1
fi


#########################


# calculate FRiP (fraction of reads in peaks)

# "ENCODE Consortium scrutinizes experiments in which the FRiP falls below 1%"
# ENCODE ATAC-seq Data Standards: ">0.3, though values greater than 0.2 are acceptable"

module add samtools/1.20

echo
echo " * samtools: $(readlink -f $(which samtools))"
echo " * samtools version: $(samtools version | grep "samtools" | head -1)"
echo

num_reads=$(samtools view -c $bam)
echo "num reads: $num_reads"

num_reads_peaks=$(samtools view -c --target-file $peaks_bed $bam)
echo "num reads in peaks: $num_reads_peaks"

frip=$(echo "(${num_reads_peaks}/${num_reads})" | bc -l | cut -c 1-4)
echo "FRiP: $frip"


#########################


# generate summary

num_peaks=$(cat "$peaks_bed" | wc -l)
echo "num peaks: $num_peaks"

# header for summary file
peaks_label="MACS3 HMMRATAC"
echo "#SAMPLE,PEAKS ${peaks_label},FRIP ${peaks_label}" > "$summary_csv"

# summarize log file
echo "${sample},${num_peaks},${frip}" >> "$summary_csv"

sleep 5

# combine all sample summaries
cat ${summary_dir}/*.${segment_name}.csv | LC_ALL=C sort -t ',' -k1,1 | uniq > "${proj_dir}/summary.${segment_name}.csv"


#########################



# end
