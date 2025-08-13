#!/bin/bash


# HMMRATAC peak calling for ATAC-seq (original Java implementation)


# script filename
script_path="${BASH_SOURCE[0]}"
script_name=$(basename "$script_path")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 4 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir sample_name BAM score \n" >&2
	exit 1
fi

# arguments
proj_dir=$1
sample=$2
bam=$3
score=$4


#########################


# settings and files

segment_name="${segment_name}-score-${score}"

summary_dir="${proj_dir}/summary"
mkdir -p "$summary_dir"
summary_csv="${summary_dir}/${sample}.${segment_name}.csv"

hmmratac_peaks_dir="${proj_dir}/peaks-HMMRATAC-score-${score}"
mkdir -p "$hmmratac_peaks_dir"
peaks_bed="${hmmratac_peaks_dir}/${sample}.peaks.bed"
summits_bed_final="${hmmratac_peaks_dir}/${sample}.summits.bed"
summits_bed_padded="${hmmratac_peaks_dir}/${sample}.summits.pad50.bed"

hmmratac_logs_dir="${proj_dir}/logs-${segment_name}"
mkdir -p "$hmmratac_logs_dir"
out_prefix="${hmmratac_logs_dir}/${sample}"
peaks_gappedpeak="${hmmratac_logs_dir}/${sample}_peaks.gappedPeak"
summits_bed_original="${hmmratac_logs_dir}/${sample}_summits.bed"

# unload all loaded modulefiles
module purge


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

if [ ! -s "${bam}.bai" ] ; then
	echo -e "\n $script_name ERROR: BAM index ${bam}.bai does not exist \n" >&2
	exit 1
fi

code_dir=$(dirname $(dirname "$script_path"))

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


# exit if output exits already

if [ -s "$peaks_bed" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 1
fi


#########################


# run HMMRATAC

hmmratac_jar="/gpfs/data/igorlab/software/HMMRATAC/HMMRATAC_V1.2.10_exe.jar"

# check that the HMMRATAC jar file is present
if [ ! -s "$hmmratac_jar" ] ; then
	echo -e "\n $script_name ERROR: jar file $hmmratac_jar does not exist \n" >&2
	exit 1
fi

# --bam <BAM> sorted BAM file containing the ATAC-seq reads
# --index <BAI> index file for the sorted BAM File
# --genome <GenomeFile> two column, tab delimited file containing genome size stats
# --blacklist <BED_File> blacklisted regions to exclude
# --threshold <double> peaks with score is >= this value will be reported (default: 30)

echo
echo " * HMMRATAC: $hmmratac_jar "
echo " * genome file: $chrom_sizes "
echo " * blacklist: $blacklist "
echo " * BAM: $bam "
echo " * out peaks gappedPeak: $peaks_gappedpeak "
echo " * out summits BED original: $summits_bed_original "
echo " * out peaks BED: $peaks_bed "
echo " * out summits BED final: $summits_bed_final "
echo

bash_cmd="
java -Xms64G -Xmx64G -jar $hmmratac_jar \
--genome $chrom_sizes \
--blacklist $blacklist \
--bam $bam \
--index ${bam}.bai \
--threshold $score \
--output $out_prefix
"
echo -e "\n CMD: $bash_cmd \n"
$bash_cmd


#########################


# check that output generated

if [ ! -s "$peaks_gappedpeak" ] ; then
	echo -e "\n $script_name ERROR: peaks $peaks_gappedpeak not generated \n" >&2
	exit 1
fi

# summits file can be empty
if [ ! -f "$summits_bed_original" ] ; then
	echo -e "\n $script_name ERROR: summits $summits_bed_original not generated \n" >&2
	exit 1
fi


#########################


# clean up

# convert peaks gappedPeak to BED format and sort
# "all skipped high coverage regions are added back to the gappedPeak output file"
# "91% of the high coverage regions that are masked in the GM12878 cells are blacklist regions"
bash_cmd="
cat $peaks_gappedpeak \
| awk '\$13 > 0' \
| cut -f 1,2,3,4,13 \
| LC_ALL=C sort -k1,1 -k2,2n \
> $peaks_bed"
echo -e "\n CMD: $bash_cmd \n"
eval "$bash_cmd"

# sort summits BED file
# "summits were likely locations of evolutionarily conserved regulatory elements"
# "could represent functional regulatory elements"
# "more robust to identify potential transcription factor binding sites"
bash_cmd="cat $summits_bed_original | LC_ALL=C sort -k1,1 -k2,2n > $summits_bed_final"
echo -e "\n CMD: $bash_cmd \n"
eval "$bash_cmd"

# generate padded summits files to call differentially accessible peaks
# https://github.com/LiuLabUB/HMMRATAC/issues/40
module add bedtools/2.27.1
bash_cmd="bedtools slop -i $summits_bed_final -g $chrom_sizes -b 50 > $summits_bed_padded"
echo -e "\n CMD: $bash_cmd \n"
eval "$bash_cmd"


#########################


# check that output generated

if [ ! -s "$peaks_gappedpeak" ] ; then
	echo -e "\n $script_name ERROR: peaks $peaks_gappedpeak not generated \n" >&2
	exit 1
fi

if [ ! -f "$summits_bed_final" ] ; then
	echo -e "\n $script_name ERROR: summits $summits_bed_final not generated \n" >&2
	exit 1
fi

if [ ! -f "$summits_bed_padded" ] ; then
	echo -e "\n $script_name ERROR: padded summits $summits_bed_padded not generated \n" >&2
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

num_peaks_unfiltered=$(cat "$peaks_gappedpeak" | wc -l)
echo "num peaks unfiltered: $num_peaks_unfiltered"

num_peaks_filtered=$(cat "$peaks_bed" | wc -l)
echo "num peaks filtered: $num_peaks_filtered"

num_summits=$(cat "$summits_bed_final" | wc -l)
echo "num summits: $num_summits"

# header for summary file
peaks_label="HMMRATAC score ${score}"
echo "#SAMPLE,PEAKS ${peaks_label},FRIP ${peaks_label}" > "$summary_csv"

# summarize log file
echo "${sample},${num_peaks_filtered},${frip}" >> "$summary_csv"


sleep 5

# combine all sample summaries
cat ${summary_dir}/*.${segment_name}.csv | LC_ALL=C sort -t ',' -k1,1 | uniq > "${proj_dir}/summary.${segment_name}.csv"


#########################



# end
