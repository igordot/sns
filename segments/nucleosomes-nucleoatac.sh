#!/bin/bash


# call nucleosomes in ATAC-Seq data using NucleoATAC


# script filename
script_name=$(basename "${BASH_SOURCE[0]}")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 4 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir sample_name threads BAM \n" >&2
	exit 1
fi

# arguments
proj_dir=$1
sample=$2
threads=$3
bam=$4


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

ref_fasta=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" REF-FASTA);

if [ ! -s "$ref_fasta" ] ; then
	echo -e "\n $script_name ERROR: FASTA $ref_fasta DOES NOT EXIST \n" >&2
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

segment_name="${segment_name}"

summary_dir="${proj_dir}/summary"
mkdir -p "$summary_dir"
summary_csv="${summary_dir}/${sample}.${segment_name}.csv"

logs_dir="${proj_dir}/logs-${segment_name}"
mkdir -p "$logs_dir"
peaks_xls="${logs_dir}/${sample}_peaks.xls"
peaks_broad="${logs_dir}/${sample}_peaks.broadPeak"
nucleoatac_basename="${logs_dir}/${sample}"
nucleoatac_occ_bg="${logs_dir}/${sample}.occ.bedgraph"
nucleoatac_nuc_smooth_bg="${logs_dir}/${sample}.nucleoatac_signal.smooth.bedgraph"

nuc_dir="${proj_dir}/nucleosomes"
mkdir -p "$nuc_dir"
peaks_bed="${nuc_dir}/${sample}.peaks.bed"
nucleoatac_occ_bed="${nuc_dir}/${sample}.occ.bed"
nucleoatac_nuc_bed="${nuc_dir}/${sample}.nucleoatac.bed"
nucleoatac_combined_bed="${nuc_dir}/${sample}.combined.bed"
nucleoatac_occ_bw="${nuc_dir}/${sample}.occ.bw"
nucleoatac_nuc_smooth_bw="${nuc_dir}/${sample}.nucleoatac.smooth.bw"


#########################


# exit if output exits already

if [ -s "$peaks_bed" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 1
fi


#########################


# MACS to call broad peaks

# MACS is part of python/2.7.3 module
module unload python
module load python/2.7.3

echo " * MACS: $(readlink -f $(which macs2)) "
echo " * MACS version: $(macs2 --version 2>&1) "
echo " * BAM: $bam "
echo " * MACS genome: $macs_genome "
echo " * MACS dir: $logs_dir "

cd "$logs_dir" || exit 1

bash_cmd="
macs2 callpeak \
--verbose 2 \
--broad --broad-cutoff 0.1 \
--format BAMPE \
--keep-dup all \
--gsize $macs_genome \
--name $sample \
--treatment $bam \
--outdir $logs_dir
"
echo -e "\n CMD: $bash_cmd \n"
$bash_cmd


#########################


# check that output generated

if [ ! -s "$peaks_xls" ] ; then
	echo -e "\n $script_name ERROR: XLS $peaks_xls NOT GENERATED \n" >&2
	exit 1
fi

if [ ! -s "$peaks_broad" ] ; then
	echo -e "\n $script_name ERROR: PEAKS $peaks_broad NOT GENERATED \n" >&2
	exit 1
fi


#########################


# generate padded broad peaks bed file

module load bedtools/2.26.0

bash_cmd="
cut -f 1,2,3,4,7 $peaks_broad \
| LC_ALL=C sort -k1,1 -k2,2n \
| bedtools slop -g $chrom_sizes -b 100 \
| bedtools merge -d 10 \
> $peaks_bed
"
echo -e "\n CMD: $bash_cmd \n"
eval "$bash_cmd"


#########################


# NucleoATAC

echo " * NucleoATAC: $(readlink -f $(which nucleoatac)) "
echo " * NucleoATAC version: $(nucleoatac --version 2>&1 | grep 'version' | grep -v 'Command') "
echo " * FASTA: $ref_fasta "
echo " * BAM: $bam "
echo " * BED: $peaks_bed "
echo " * OUT BASE: $nucleoatac_basename "

bash_cmd="
nucleoatac run \
--cores $threads \
--fasta $ref_fasta \
--bam $bam \
--bed $peaks_bed \
--out $nucleoatac_basename
"
echo -e "\n CMD: $bash_cmd \n"
$bash_cmd


#########################


# check that output generated

if [ ! -s "${nucleoatac_basename}.nucpos.bed.gz" ] ; then
	echo -e "\n $script_name ERROR: BED ${nucleoatac_basename}.nucpos.bed.gz NOT GENERATED \n" >&2
	exit 1
fi

if [ ! -s "${nucleoatac_basename}.nucmap_combined.bed.gz" ] ; then
	echo -e "\n $script_name ERROR: BED ${nucleoatac_basename}.nucmap_combined.bed.gz NOT GENERATED \n" >&2
	exit 1
fi

gunzip ${nucleoatac_occ_bg}.gz

if [ ! -s "$nucleoatac_occ_bg" ] ; then
	echo -e "\n $script_name ERROR: BEDGRAPH $nucleoatac_occ_bg NOT GENERATED \n" >&2
	exit 1
fi

gunzip ${nucleoatac_nuc_smooth_bg}.gz

if [ ! -s "$nucleoatac_nuc_smooth_bg" ] ; then
	echo -e "\n $script_name ERROR: BEDGRAPH $nucleoatac_nuc_smooth_bg NOT GENERATED \n" >&2
	exit 1
fi


#########################


# convert bedGraphs to bigWigs for easier consumption

module load kentutils/329

echo " * OCC bedGraph: $nucleoatac_occ_bg "
echo " * OCC bigWig: $nucleoatac_occ_bw "

bw_cmd="bedGraphToBigWig $nucleoatac_occ_bg $chrom_sizes $nucleoatac_occ_bw"
echo -e "\n CMD: $bw_cmd \n"
$bw_cmd

echo " * NUC bedGraph: $nucleoatac_nuc_smooth_bg "
echo " * NUC bigWig: $nucleoatac_nuc_smooth_bw "

bw_cmd="bedGraphToBigWig $nucleoatac_nuc_smooth_bg $chrom_sizes $nucleoatac_nuc_smooth_bw"
echo -e "\n CMD: $bw_cmd \n"
$bw_cmd


#########################


# move or delete files

# EPS plots
mv -v ${nucleoatac_basename}.*.eps ${nuc_dir}/

# converting to png results in quality loss
# convert -density 300 -flatten $f ${f%.*}.png;
# mogrify -format png *.eps

gunzip ${nucleoatac_basename}.*.bed.gz

# positions of nucleosomes estimated from the *.occ.bedgraph.gz (low resolution)
mv -v "${nucleoatac_basename}.occpeaks.bed" "$nucleoatac_occ_bed"
# positions of nucleosomes estimated from the *.nucleoatac_signal.bedgraph.gz track (higher resolution)
mv -v "${nucleoatac_basename}.nucpos.bed" "$nucleoatac_nuc_bed"
# merge of the other two, with the positions from the nucleoatac_signal favored when there is an overlap
mv -v "${nucleoatac_basename}.nucmap_combined.bed" "$nucleoatac_combined_bed"

# delete files that are not needed
rm -rfv ${nucleoatac_basename}.*.bed.gz.tbi
rm -rfv ${nucleoatac_basename}.occ.lower_bound.bedgraph*
rm -rfv ${nucleoatac_basename}.occ.upper_bound.bedgraph*
rm -rfv ${nucleoatac_basename}.ins.bedgraph*

rm -fv ${nucleoatac_occ_bg}.gz.tbi
rm -fv ${nucleoatac_nuc_smooth_bg}.gz.tbi


#########################


# generate summary

nucleosomes_occ_num=$(cat "$nucleoatac_occ_bed" | wc -l)
echo "nucleosomes occ: $nucleosomes_occ_num"

nucleosomes_nuc_num=$(cat "$nucleoatac_nuc_bed" | wc -l)
echo "nucleosomes nuc: $nucleosomes_nuc_num"

nucleosomes_combined_num=$(cat "$nucleoatac_combined_bed" | wc -l)
echo "nucleosomes combined: $nucleosomes_combined_num"

# header for summary file
echo "#SAMPLE,nucleosomes occ,nucleosomes nuc,nucleosomes combined" > "$summary_csv"

# summarize log file
echo "${sample},${nucleosomes_occ_num},${nucleosomes_nuc_num},${nucleosomes_combined_num}" >> "$summary_csv"

sleep 30

# combine all sample summaries
cat ${summary_dir}/*.${segment_name}.csv | LC_ALL=C sort -t ',' -k1,1 | uniq > "${proj_dir}/summary.${segment_name}.csv"


#########################



# end
