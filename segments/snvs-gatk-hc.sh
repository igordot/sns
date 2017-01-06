#!/bin/bash


# call variants with GATK HaplotypeCaller


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

vcf_dir="${proj_dir}/VCF-GATK-HC"
mkdir -p "$vcf_dir"
vcf_original="${vcf_dir}/${sample}.original.vcf"
vcf_fixed="${vcf_dir}/${sample}.vcf"


#########################


# exit if output exits already

if [ -s "$vcf_original" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 1
fi

if [ -s "$vcf_fixed" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 1
fi


#########################


# GATK settings

module unload java
module load java/1.8

# command
gatk_jar="/ifs/home/id460/software/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar"
gatk_cmd="java -Xms16G -Xmx16G -jar ${gatk_jar}"

# error log (DEBUG, INFO (default), WARN, ERROR, FATAL, OFF)
gatk_log_level_arg="--logging_level ERROR"

if [ ! -s "$gatk_jar" ] ; then
	echo -e "\n $script_name ERROR: GATK $gatk_jar DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# GATK HaplotypeCaller

echo " * GATK: $(readlink -f $gatk_jar) "
echo " * GATK version: $($gatk_cmd --version) "
echo " * BAM: $bam "
echo " * INTERVALS: $bed "
echo " * VCF original: $vcf_original "
echo " * VCF fixed: $vcf_fixed "

gatk_hc_cmd="
$gatk_cmd -T HaplotypeCaller -dt NONE $gatk_log_level_arg \
-nct $threads \
--max_alternate_alleles 3 \
--standard_min_confidence_threshold_for_calling 50 \
--reference_sequence $ref_fasta \
--intervals $bed \
--interval_padding 10 \
--input_file $bam \
--out $vcf_original
"
echo -e "\n CMD: $gatk_hc_cmd \n"
$gatk_hc_cmd


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

# 1) GATK HC is defining the AD field as "Number=." (VCF 4.1 specification) rather than "Number=R" (VCF 4.2 specification)
# CMD: sed 's/AD,Number=./AD,Number=R/g' $vcf_original
# fixed in GATK 3.6

# 1) split multi-allelic variants calls into separate lines (uses VCF 4.2 specification)
# 2) perform indel left-normalization (start position of a variant should be shifted to the left until it is no longer possible to do so)
# 3) depth filter

fix_vcf_cmd="
cat $vcf_original \
| bcftools norm --multiallelics -both --output-type v - \
| bcftools norm --fasta-ref $ref_fasta --output-type v - \
| bcftools view --exclude 'DP<5' --output-type v > $vcf_fixed
"
echo -e "\n CMD: $fix_vcf_cmd \n"
eval "$fix_vcf_cmd"


#########################


# check that output generated

if [ ! -s "$vcf_fixed" ] ; then
	echo -e "\n $script_name ERROR: VCF $vcf_fixed NOT GENERATED \n" >&2
	exit 1
fi


#########################


# annotate

annot_cmd="bash ${code_dir}/segments/annot-annovar.sh $proj_dir $sample $vcf_fixed"
echo -e "\n CMD: $annot_cmd \n"
($annot_cmd)


#########################



# end
