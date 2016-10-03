#!/bin/bash


# Split reads that contain Ns in their CIGAR string with GATK SplitNCigarReads


# script filename
script_name=$(basename "${BASH_SOURCE[0]}")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 3 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir sample_name BAM \n" >&2
	exit 1
fi

# arguments
proj_dir=$1
sample=$2
bam=$3


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

ref_fasta=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" REF-FASTA);

if [ ! -s "$ref_fasta" ] ; then
	echo -e "\n $script_name ERROR: FASTA $ref_fasta DOES NOT EXIST \n" >&2
	exit 1
fi

ref_dict=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" REF-DICT);

if [ ! -s "$ref_dict" ] ; then
	echo -e "\n $script_name ERROR: DICT $ref_dict DOES NOT EXIST \n" >&2
	exit 1
fi

gtf=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" REF-GTF);

if [ ! -s "$gtf" ] || [ ! "$gtf" ] ; then
	echo -e "\n $script_name ERROR: GTF $gtf DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# settings and files

samples_csv="${proj_dir}/samples.${segment_name}.csv"

# account for both dedup (.dd.bam) non-dedup (.bam) input BAMs
bam_base=$(basename "$bam")
bam_base=${bam_base/%.bam/}

bam_split_dir="${proj_dir}/BAM-SPLIT"
mkdir -p "$bam_split_dir"
bam_split="${bam_split_dir}/${bam_base}.bam"


#########################


# exit if output exits already

if [ -s "$bam_split" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 1
fi


#########################


# GATK settings

module unload java
module load java/1.8
module load r/3.3.0

# command
gatk_jar="/ifs/home/id460/software/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar"
gatk_cmd="java -Xms16G -Xmx16G -jar ${gatk_jar}"

if [ ! -s "$gatk_jar" ] ; then
	echo -e "\n $script_name ERROR: GATK $gatk_jar DOES NOT EXIST \n" >&2
	exit 1
fi

# error log (blank for troubleshooting)
gatk_log_level_arg="--logging_level ERROR"


#########################


# GATK SplitNCigarReads

echo " * GATK: $(readlink -f $gatk_jar) "
echo " * GATK version: $($gatk_cmd --version) "
echo " * BAM in: $bam "
echo " * BAM out: $bam_split "

gatk_split_cmd="
$gatk_cmd -T SplitNCigarReads $gatk_log_level_arg \
--unsafe ALLOW_N_CIGAR_READS
--reference_sequence $ref_fasta \
--input_file $bam \
--out $bam_split
"
echo "CMD: $gatk_split_cmd"
$gatk_split_cmd


#########################


# check that output generated

if [ ! -s "$bam_split" ] ; then
	echo -e "\n $script_name ERROR: BAM $bam_split NOT GENERATED \n" >&2
	exit 1
fi


#########################


# exons BED file (needed for future steps)

found_bed=$(find $proj_dir -maxdepth 1 -type f -name "*.bed" | head -1)
bed=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" EXP-BED $found_bed);

# generate exons BED file if doesn't exist already
if [ ! -s "$bed" ] ; then

	module load bedtools/2.26.0

	cat "$gtf" \
	| awk -F $'\t' '$3 == "exon" && $5 > $4' \
	| cut -f 1,4,5 \
	| LC_ALL=C sort -u -k1,1 -k2,2n \
	| bedtools merge \
	> "${proj_dir}/exons.bed"

	found_bed=$(find $proj_dir -maxdepth 1 -type f -name "*.bed" | head -1)
	bed=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" EXP-BED $found_bed);

fi


#########################


# add sample and BAM to sample sheet
echo "${sample},${bam_split}" >> "$samples_csv"

sleep 30

# sort and remove duplicates in place in sample sheet
LC_ALL=C sort -t ',' -k1,1 -u -o "$samples_csv" "$samples_csv"

sleep 30


#########################



# end
