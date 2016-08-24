#!/bin/bash


# bedtools generate BigWig from BAM


# script filename
script_name=$(basename "${BASH_SOURCE[0]}")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 3 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name [project dir] [sample] [BAM] \n" >&2
	exit 1
fi

# arguments
PROJ_DIR=$1
SAMPLE=$2
BAM=$3


#########################


# check that inputs exist

if [ ! -d "$PROJ_DIR" ] || [ ! "$PROJ_DIR" ] ; then
	echo -e "\n $script_name ERROR: DIR $PROJ_DIR DOES NOT EXIST \n" >&2
	exit 1
fi

if [ ! -s "$BAM" ] || [ ! "$BAM" ] ; then
	echo -e "\n $script_name ERROR: BAM $BAM DOES NOT EXIST \n" >&2
	exit 1
fi

CODE_DIR=$(dirname "$(dirname "${BASH_SOURCE[0]}")")
REF_CHROMSIZES=$(bash ${CODE_DIR}/scripts/get-set-setting.sh "${PROJ_DIR}/settings.txt" REF-CHROMSIZES);

if [ ! -s "$REF_CHROMSIZES" ] ; then
	echo -e "\n $script_name ERROR: CHROMSIZES $REF_CHROMSIZES DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# output

BIGWIG_DIR="${PROJ_DIR}/BIGWIG"
mkdir -p "$BIGWIG_DIR"
BEDGRAPH="${BIGWIG_DIR}/${SAMPLE}.norm.bedgraph"
BIGWIG="${BIGWIG_DIR}/${SAMPLE}.norm.bw"


#########################


# exit if output exits already

if [ -s "$BIGWIG" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $SAMPLE \n" >&2
	exit 1
fi


#########################


# BAM to BedGraph with bedtools

module unload gcc
module load samtools/1.3
module load bedtools/2.22.0
module load kentutils/329

# determine normalization scaling factor
MAPPED_READS=$(samtools view -c $BAM)
SCALE=`echo "10000000/${MAPPED_READS}" | bc -l`
SCALE=${SCALE:0:8}

echo " * samtools: $(readlink -f $(which samtools)) "
echo " * bedtools: $(readlink -f $(which bedtools)) "
echo " * bedGraphToBigWig: $(readlink -f $(which bedGraphToBigWig)) "
echo " * CHROM SIZES: $REF_CHROMSIZES "
echo " * BAM: $BAM "
echo " * MAPPED READS: $MAPPED_READS "
echo " * SCALING FACTOR: $SCALE "
echo " * BEDGRAPH: $BEDGRAPH "
echo " * BIGWIG: $BIGWIG "

CMD="bedtools genomecov -split -bg -g $REF_CHROMSIZES -ibam $BAM -scale $SCALE > $BEDGRAPH"
echo "CMD: $CMD"
eval $CMD


#########################


# check that output generated

if [ ! -s "$BEDGRAPH" ] ; then
	echo -e "\n $script_name ERROR: BEDGRAPH $BEDGRAPH NOT GENERATED \n" >&2
	exit 1
fi


#########################


# BadGraph to BigWig with UCSC tools

CMD="bedGraphToBigWig $BEDGRAPH $REF_CHROMSIZES $BIGWIG"
echo "CMD: $CMD"
eval $CMD

# delete bedgraph
rm -f $BEDGRAPH


#########################


# check that output generated

if [ ! -s "$BIGWIG" ] ; then
	echo -e "\n $script_name ERROR: BIGWIG $BIGWIG NOT GENERATED \n" >&2
	exit 1
fi


#########################



# end
