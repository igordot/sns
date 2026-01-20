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
	echo -e "\n USAGE: $script_name project_dir peak_type q_value sample_name treatment_bam [control_bam] \n" >&2
	if [ $# -gt 0 ] ; then echo -e "\n ARGS: $* \n" >&2 ; fi
	exit 1
fi

# arguments
proj_dir=$(readlink -f "$1")
peak_type=$2
q_value=$3
sample=$4
bam_treat=$5
bam_control=$6


#########################


# settings and files

# adjust segment label based on run type and presence of control sample
if [ "$peak_type" == "atac" ] ; then
	settings_label="${peak_type}"
elif [ "$peak_type" == "narrow" ] || [ "$peak_type" == "broad" ] ; then
	settings_label="${peak_type}"
	if [ -z "$bam_control" ] ; then
		settings_label="${settings_label}-nocontrol"
	fi
else
	echo -e "\n $script_name ERROR: unknown peak type $peak_type \n" >&2
	exit 1
fi

# add q-value to segment label
settings_label="${settings_label}-q-${q_value}"

# add settings to the segment label
segment_name="${segment_name}-${settings_label}"

summary_dir="${proj_dir}/summary"
mkdir -p "$summary_dir"
summary_csv="${summary_dir}/${sample}.${segment_name}.csv"

macs_peaks_dir="${proj_dir}/peaks-MACS3-${settings_label}"
mkdir -p "$macs_peaks_dir"
peaks_bed="${macs_peaks_dir}/${sample}.peaks.bed"
summits_bed="${macs_peaks_dir}/${sample}.summits.bed"

macs_logs_dir="${proj_dir}/logs-${segment_name}"
mkdir -p "$macs_logs_dir"
peaks_xls="${macs_logs_dir}/${sample}_peaks.xls"
peaks_file="${macs_logs_dir}/${sample}_peaks.narrowPeak"
summits_file="${macs_logs_dir}/${sample}_summits.bed"
model_r="${macs_logs_dir}/${sample}_model.r"
macs_bdg_treat="${macs_logs_dir}/${sample}_treat_pileup.bdg"
macs_bdg_treat_sort="${macs_logs_dir}/${sample}_treat_pileup.sort.bdg"
macs_bdg_control="${macs_logs_dir}/${sample}_control_lambda.bdg"

bigwig_dir="${proj_dir}/BIGWIG"
mkdir -p "$bigwig_dir"
macs_bw="${bigwig_dir}/${sample}.macs3.bw"

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

if [ ! -s "$bam_treat" ] ; then
	echo -e "\n $script_name ERROR: BAM $bam_treat does not exist \n" >&2
	exit 1
fi

code_dir=$(dirname $(dirname "$script_path"))

genome_dir=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" GENOME-DIR);

if [ ! -d "$genome_dir" ] ; then
	echo -e "\n $script_name ERROR: genome dir $genome_dir does not exist \n" >&2
	exit 1
fi

chrom_sizes=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" REF-CHROMSIZES);

if [ ! -s "$chrom_sizes" ] ; then
	echo -e "\n $script_name ERROR: chrom sizes $chrom_sizes does not exist \n" >&2
	exit 1
fi

blacklist=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" REF-BLACKLIST);

if [ ! -s "$blacklist" ] ; then
	echo -e "\n $script_name ERROR: blacklist $blacklist does not exist \n" >&2
	exit 1
fi


#########################


# MACS parameters

# MACS-style genome abbreviation (keep first two characters of build name)
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
	analysis_params="--nomodel --shift -100 --extsize 200 --qvalue $q_value"
elif [ "$peak_type" == "narrow" ] ; then
	analysis_params="--qvalue $q_value"
elif [ "$peak_type" == "broad" ] ; then
	analysis_params="--broad --broad-cutoff $q_value"
	peaks_file="${macs_logs_dir}/${sample}_peaks.broadPeak"
else
	echo
fi


#########################


# MACS

# MACS3 is part of condaenvs/2023/macs3 module
module add condaenvs/2023/macs3

# check that the binary can be run
if ! macs3 --version >/dev/null 2>&1; then
	echo -e "\n $script_name ERROR: macs3 cannot be executed \n" >&2
	exit 1
fi
if [ -z "$(macs3 --version 2>/dev/null)" ]; then
	echo -e "\n $script_name ERROR: macs3 cannot be executed \n" >&2
	exit 1
fi

echo
echo " * MACS: $(readlink -f $(which macs3)) "
echo " * MACS version: $(macs3 --version 2>&1) "
echo " * Python: $(readlink -f $(which python)) "
echo " * Python version: $(python --version 2>&1) "
echo " * BAM treatment: $bam_treat "
echo " * BAM control: $bam_control "
echo " * MACS genome: $macs_genome "
echo " * MACS logs dir: $macs_logs_dir "
echo " * MACS peaks dir: $macs_peaks_dir "
echo " * peak type: $peak_type "
echo " * peak cutoff: $q_value "
echo

cd "$macs_logs_dir" || exit 1

bash_cmd="
macs3 callpeak \
--verbose 2 \
--bdg --SPMR --call-summits \
--keep-dup all \
$analysis_params \
--gsize $macs_genome \
--name $sample \
$bam_param \
--outdir $macs_logs_dir
"
echo -e "\n CMD: $bash_cmd \n"
$bash_cmd


#########################


# check that output generated

if [ ! -s "$peaks_file" ] ; then
	echo -e "\n $script_name ERROR: peaks $peaks_file not generated \n" >&2
	exit 1
fi

if [ ! -s "$peaks_xls" ] ; then
	echo -e "\n $script_name ERROR: XLS $peaks_xls not generated \n" >&2
	exit 1
fi


#########################


# generate an image about the model based on the data
# --nomodel will bypass building the shifting model

module add r/4.1.2

if [ -s "$model_r" ] ; then

	echo
	echo " * R: $(readlink -f $(which R)) "
	echo " * R version: $(R --version | head -1) "
	echo " * Rscript: $(readlink -f $(which Rscript)) "
	echo " * Rscript version: $(Rscript --version 2>&1) "
	echo

	bash_cmd="Rscript $model_r"
	echo -e "\n CMD: $bash_cmd \n"
	$bash_cmd

fi

sleep 5


#########################


# generate a blacklist-filtered BED file

module purge
module add bedtools/2.30.0

# check that the binary can be run
if ! bedtools --version >/dev/null 2>&1; then
	echo -e "\n $script_name ERROR: bedtools cannot be executed \n" >&2
	exit 1
fi

echo
echo " * bedtools: $(readlink -f $(which bedtools)) "
echo " * bedtools version: $(bedtools --version) "
echo " * blacklist: $blacklist "
echo

# keep only peaks that do not overlap blacklist regions
bash_cmd="
cut -f 1,2,3,4,7 $peaks_file \
| bedtools intersect -v -a stdin -b $blacklist \
| LC_ALL=C sort -k1,1 -k2,2n \
> $peaks_bed
"
echo -e "\n CMD: $bash_cmd \n"
eval "$bash_cmd"

# keep only summits of filtered peaks (generated for narrow peaks)
if [ -s "$summits_file" ] ; then
	bash_cmd="
	bedtools intersect -wa -a $summits_file -b $peaks_bed \
	| LC_ALL=C sort -k1,1 -k2,2n \
	| uniq \
	> $summits_bed
	"
	echo -e "\n CMD: $bash_cmd \n"
	eval "$bash_cmd"
fi


#########################


# check that output generated

if [ ! -s "$peaks_bed" ] ; then
	echo -e "\n $script_name ERROR: peaks $peaks_bed not generated \n" >&2
	exit 1
fi


#########################


# generate bigWig file if not already generated

if [ ! -s "$macs_bdg_treat" ] ; then
	echo -e "\n $script_name ERROR: bedGraph $macs_bdg_treat not generated \n" >&2
	exit 1
fi

# bedGraphToBigWig is part of UCSC Genome Browser Group's suite of programs
# ucscutils/374 requires mariadb/5.5.64 to be loaded
module add ucscutils/398

# check that the binary can be run
if ! bedGraphToBigWig 2>&1 | grep -q "usage"; then
	echo -e "\n $script_name ERROR: bedGraphToBigWig cannot be executed \n" >&2
	exit 1
fi

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


# check that output generated

if [ ! -s "$macs_bw" ] ; then
	echo -e "\n $script_name ERROR: bigWig $macs_bw not generated \n" >&2
	exit 1
fi


#########################


# calculate FRiP (fraction of reads in peaks)

# "ENCODE Consortium scrutinizes experiments in which the FRiP falls below 1%"
# ENCODE ATAC-seq Data Standards: ">0.3, though values greater than 0.2 are acceptable"

module add samtools/1.20

# check that the binary can be run
if ! samtools --version >/dev/null 2>&1; then
	echo -e "\n $script_name ERROR: samtools cannot be executed \n" >&2
	exit 1
fi

echo
echo " * samtools: $(readlink -f $(which samtools))"
echo " * samtools version: $(samtools version | grep "samtools" | head -1)"
echo

num_reads=$(samtools view -c $bam_treat)
echo "num reads: $num_reads"

num_reads_peaks=$(samtools view -c --target-file $peaks_bed $bam_treat)
echo "num reads in peaks: $num_reads_peaks"

frip=$(echo "(${num_reads_peaks}/${num_reads})" | bc -l | cut -c 1-4)
echo "FRiP: $frip"


#########################


# generate summary

num_peaks_unfiltered=$(cat "$peaks_file" | cut -f 1-3 | uniq | wc -l)
echo "num peaks unfiltered: $num_peaks_unfiltered"

num_peaks_filtered=$(cat "$peaks_bed" | cut -f 1-3 | uniq | wc -l)
echo "num peaks filtered: $num_peaks_filtered"

num_subpeaks=$(cat "$peaks_bed" | wc -l)
echo "num subpeaks: $num_subpeaks"

num_summits=$(cat "$summits_bed" | wc -l)
echo "num summits: $num_summits"

# header for summary file
peaks_label="MACS3 ${peak_type} q ${q_value}"
echo "#SAMPLE,PEAKS ${peaks_label},FRIP ${peaks_label}" > "$summary_csv"

# summarize log file
echo "${sample},${num_peaks_filtered},${frip}" >> "$summary_csv"

sleep 5

# combine all sample summaries
cat ${summary_dir}/*.${segment_name}.csv | LC_ALL=C sort -t ',' -k1,1 | uniq > "${proj_dir}/summary.${segment_name}.csv"


#########################



# end
