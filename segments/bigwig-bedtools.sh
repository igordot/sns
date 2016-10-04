#!/bin/bash


# bedtools generate bigWig from BAM


# script filename
script_name=$(basename "${BASH_SOURCE[0]}")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 3 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir sample BAM \n" >&2
	exit 1
fi

# arguments
proj_dir=$1
sample=$2
bam=$3


#########################


# check that inputs exist

if [ ! -d "$proj_dir" ] || [ ! "$proj_dir" ] ; then
	echo -e "\n $script_name ERROR: DIR $proj_dir DOES NOT EXIST \n" >&2
	exit 1
fi

if [ ! -s "$bam" ] || [ ! "$bam" ] ; then
	echo -e "\n $script_name ERROR: BAM $bam DOES NOT EXIST \n" >&2
	exit 1
fi

code_dir=$(dirname "$(dirname "${BASH_SOURCE[0]}")")
chrom_sizes=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" REF-CHROMSIZES);

if [ ! -s "$chrom_sizes" ] ; then
	echo -e "\n $script_name ERROR: CHROM SIZES $chrom_sizes DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# output

BIGWIG_DIR="${proj_dir}/BIGWIG"
mkdir -p "$BIGWIG_DIR"
bedgraph="${BIGWIG_DIR}/${sample}.norm.bedgraph"
bigwig="${BIGWIG_DIR}/${sample}.norm.bw"


#########################


# exit if output exits already

if [ -s "$bigwig" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 1
fi


#########################


# bam to BedGraph with bedtools

module load samtools/1.3
module load bedtools/2.26.0
module load kentutils/329

# determine normalization scaling factor
mapped_reads=$(samtools view -c $bam)
scale=$(echo "10000000/${mapped_reads}" | bc -l)
scale=${scale:0:8}

echo " * samtools: $(readlink -f $(which samtools)) "
echo " * bedtools: $(readlink -f $(which bedtools)) "
echo " * bedGraphToBigWig: $(readlink -f $(which bedGraphToBigWig)) "
echo " * chrom sizes: $chrom_sizes "
echo " * BAM: $bam "
echo " * mapped reads: $mapped_reads "
echo " * scaling factor: $scale "
echo " * BedGraph: $bedgraph "
echo " * bigWig: $bigwig "

# bedGraphToBigWig error:
# ".bedgraph is not case-sensitive sorted at line 2771141.  Please use "sort -k1,1 -k2,2n" with LC_COLLATE=C,  or bedSort and try again."

bg_cmd="bedtools genomecov -split -bg -g $chrom_sizes -ibam $bam -scale $scale | LC_ALL=C sort -k1,1 -k2,2n > $bedgraph"
echo -e "\n CMD: $bg_cmd \n"
eval "$bg_cmd"


#########################


# check that output generated

if [ ! -s "$bedgraph" ] ; then
	echo -e "\n $script_name ERROR: BEDGRAPH $bedgraph NOT GENERATED \n" >&2
	exit 1
fi


#########################


# BadGraph to BigWig with UCSC tools

bw_cmd="bedGraphToBigWig $bedgraph $chrom_sizes $bigwig"
echo -e "\n CMD: $bw_cmd \n"
eval "$bw_cmd"

# delete bedgraph
rm -fv $bedgraph


#########################


# check that output generated

if [ ! -s "$bigwig" ] ; then
	echo -e "\n $script_name ERROR: BIGWIG $bigwig NOT GENERATED \n" >&2
	exit 1
fi


#########################



# end
