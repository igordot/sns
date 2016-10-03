#!/bin/bash


# call variants with LoFreq


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

bed=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" EXP-BED $found_bed);

if [ ! -s "$bed" ] ; then
	echo -e "\n $script_name ERROR: BED $bed DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# settings and files

# summary_dir="${proj_dir}/summary"
# mkdir -p "$summary_dir"
# summary_csv="${summary_dir}/${sample}.${segment_name}.csv"

vcf_dir="${proj_dir}/VCF-LoFreq"
mkdir -p "$vcf_dir"
vcf_original="${vcf_dir}/${sample}.original.vcf"
vcf_fixed="${vcf_dir}/${sample}.vcf"


#########################


# exit if output exits already

if [ -s "$vcf_fixed" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 1
fi


#########################


# create padded bed if needed

bed_padded="${bed}.pad10"

if [ ! -s $bed_padded ] ; then
	module load bedtools/2.26.0
	bedtools slop -i $bed -g $chrom_sizes -b 10 > $bed_padded
	sleep 30
fi

if [ ! -s $bed_padded ] ; then
	echo -e "\n $script_name ERROR: BED $bed_padded DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# LoFreq

lofreq_bin="/ifs/home/id460/bin/lofreq"

echo " * LoFreq: $(readlink -f $(which $lofreq_bin)) "
echo " * LoFreq version: $($lofreq_bin version 2>&1 | head -1) "
echo " * BAM: $bam "
echo " * BED: $bed_padded "
echo " * VCF original: $vcf_original "
echo " * VCF fixed: $vcf_fixed "

lofreq_cmd="
$lofreq_bin call-parallel --call-indels \
--pp-threads $threads \
--ref $ref_fasta \
--bed $bed_padded \
--out $vcf_original \
$bam
"
echo "CMD: $lofreq_cmd"
$lofreq_cmd


#########################


# check that output generated

if [ ! -s "$vcf_original" ] ; then
	echo -e "\n $script_name ERROR: VCF $vcf_original NOT GENERATED \n" >&2
	exit 1
fi


#########################


# adjust the vcf for annovar compatibility (http://www.openbioinformatics.org/annovar/annovar_vcf.html)

module unload samtools
module load samtools/1.3

# create indexed VCF file
# warnings: "contig 'chr1' is not defined in the header. (Quick workaround: index the file with tabix.)"

# bgzip VCF file so it can be indexed
bgzip_cmd="bgzip -c $vcf_original > ${vcf_original}.bgz"
echo "CMD: $bgzip_cmd"
eval "$bgzip_cmd"

# index VCF file
bcf_index_cmd="bcftools index ${vcf_original}.bgz"
echo "CMD: $bcf_index_cmd"
eval "$bcf_index_cmd"

# 1) split multi-allelic variants calls into separate lines (uses VCF 4.2 specification)
# 2) perform indel left-normalization (start position of a variant should be shifted to the left until it is no longer possible to do so)
# 3) depth filter

fix_vcf_cmd="
bcftools norm --multiallelics -both --output-type v ${vcf_original}.bgz \
| bcftools norm --fasta-ref $ref_fasta --output-type v - \
| bcftools view --exclude 'DP<5' --output-type v > $vcf_fixed
"
echo "CMD: $fix_vcf_cmd"
eval "$fix_vcf_cmd"

# add FORMAT and sample ID columns so VCF is compatible with some other scripts
sed -i "s/FILTER\tINFO/FILTER\tINFO\tFORMAT\t${sample}/g" $vcf_fixed


#########################


# check that output generated

if [ ! -s "$vcf_fixed" ] ; then
	echo -e "\n $script_name ERROR: VCF $vcf_fixed NOT GENERATED \n" >&2
	exit 1
fi


#########################



# end
