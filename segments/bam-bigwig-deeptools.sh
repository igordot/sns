#!/bin/bash


# deepTools generate BigWig from BAM


# script filename
script_name=$(basename "${BASH_SOURCE[0]}")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 4 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name [project dir] [sample] [threads] [BAM] \n" >&2
	exit 1
fi

# arguments
PROJ_DIR=$1
SAMPLE=$2
THREADS=$3
BAM=$4


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


#########################


# output

BIGWIG_DIR="${PROJ_DIR}/BIGWIG"
mkdir -p "$BIGWIG_DIR"
BIGWIG="${BIGWIG_DIR}/${SAMPLE}.bin1.rpkm.bw"

# BIGWIG_RPKM_BIN1=${BIGWIG_RPKM_DIR}/${ID}.bin1.rpkm.bw
# BIGWIG_RPKM_BIN10=${BIGWIG_RPKM_DIR}/${ID}.bin10.smooth50.rpkm.bw


#########################


# exit if output exits already

if [ -s "$BIGWIG" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $SAMPLE \n" >&2
	exit 1
fi


#########################


# generate bigwig using deeptools

module unload python
module load python/2.7.3
module load deeptools/2.2.4
module load kentutils/329

echo " * bamCoverage: $(readlink -f $(which bamCoverage)) "
echo " * BAM: $BAM "
echo " * BIGWIG: $BIGWIG "

CMD="
bamCoverage \
--verbose \
--numberOfProcessors $THREADS \
--binSize 1 \
--normalizeUsingRPKM \
--outFileFormat bigwig \
--bam $BAM \
--outFileName $BIGWIG
"
echo "CMD: $CMD"
eval "$CMD"

# echo " * BIGWIG: $BIGWIG_RPKM_BIN10 "

# CMD="
# bamCoverage \
# --verbose \
# --numberOfProcessors $THREADS \
# --binSize 10 \
# --fragmentLength 1 \
# --smoothLength 50 \
# --normalizeUsingRPKM \
# --outFileFormat bigwig \
# --bam $BAM \
# --outFileName $BIGWIG_RPKM_BIN10
# "
# echo "CMD: $CMD"
# eval "$CMD"


#########################


# check that output generated

if [ ! -s "$BIGWIG" ] ; then
	echo -e "\n $script_name ERROR: BIGWIG $BIGWIG NOT GENERATED \n" >&2
	exit 1
fi


#########################



# end
