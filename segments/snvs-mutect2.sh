#!/bin/bash


# call variants with MuTect2 (part of GATK)


# script filename
script_name=$(basename "${BASH_SOURCE[0]}")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 5 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir tumor_sample_name tumor_bam normal_sample_name normal_bam \n" >&2
	exit 1
fi

# arguments
proj_dir=$1
sample_t=$2
bam_t=$3
sample_n=$4
bam_n=$5


#########################


# check that inputs exist

if [ ! -d "$proj_dir" ] ; then
	echo -e "\n $script_name ERROR: PROJ DIR $proj_dir DOES NOT EXIST \n" >&2
	exit 1
fi

if [ ! -s "$bam_t" ] ; then
	echo -e "\n $script_name ERROR: BAM $bam_t DOES NOT EXIST \n" >&2
	exit 1
fi

if [ ! -s "$bam_n" ] ; then
	echo -e "\n $script_name ERROR: BAM $bam_n DOES NOT EXIST \n" >&2
	exit 1
fi

code_dir=$(dirname "$(dirname "${BASH_SOURCE[0]}")")

genome_dir=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" GENOME-DIR)

if [ ! -d "$genome_dir" ] ; then
	echo -e "\n $script_name ERROR: GENOME DIR $genome_dir DOES NOT EXIST \n" >&2
	exit 1
fi

ref_fasta=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" REF-FASTA)

if [ ! -s "$ref_fasta" ] ; then
	echo -e "\n $script_name ERROR: FASTA $ref_fasta DOES NOT EXIST \n" >&2
	exit 1
fi

bed=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" EXP-TARGETS-BED)

if [ ! -s "$bed" ] ; then
	echo -e "\n $script_name ERROR: BED $bed DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# settings and files

sample="${sample_t}:${sample_n}"

vcf_dir="${proj_dir}/VCF-MuTect2"
mkdir -p "$vcf_dir"
vcf_original="${vcf_dir}/${sample_t}-${sample_n}.original.vcf"
idx_original="${vcf_original}.idx"
vcf_fixed="${vcf_dir}/${sample_t}-${sample_n}.vcf"


#########################


# skip to annotation if output exists already

annot_cmd="bash ${code_dir}/segments/annot-annovar.sh $proj_dir $sample $vcf_fixed"

if [ -s "$vcf_fixed" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	echo -e "\n CMD: $annot_cmd \n"
	($annot_cmd)
	exit 1
fi

if [ -s "$vcf_original" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 1
fi


#########################


# GATK settings

module unload java
module load java/1.8

# command
gatk_jar="/ifs/home/id460/software/GenomeAnalysisTK/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar"
gatk_cmd="java -Xms16G -Xmx16G -jar ${gatk_jar}"

if [ ! -s "$gatk_jar" ] ; then
	echo -e "\n $script_name ERROR: GATK $gatk_jar DOES NOT EXIST \n" >&2
	exit 1
fi

# error log (DEBUG, INFO (default), WARN, ERROR, FATAL, OFF)
gatk_log_level_arg="--logging_level WARN"

# known variants (may vary greatly for each genome)
genome_build=$(basename "$genome_dir")
if [[ "$genome_build" == "hg19" ]] ; then
	cosmic_vcf="${genome_dir}/CosmicCodingMuts_v73.hg19.vcf"
	dbsnp_vcf="${genome_dir}/gatk-bundle/dbsnp_138.hg19.vcf"
	mt_var_arg="--dbsnp $dbsnp_vcf --cosmic $cosmic_vcf"
elif [[ "$genome_build" == "mm10" ]] ; then
	dbsnp_vcf="${genome_dir}/dbSNP/dbsnp.146.vcf"
	mt_var_arg="--dbsnp $dbsnp_vcf"
else
	mt_var_arg=""
fi


#########################


# GATK MuTect2

echo " * GATK: $(readlink -f $gatk_jar) "
echo " * GATK version: $($gatk_cmd --version) "
echo " * sample T : $sample_t "
echo " * BAM T : $bam_t "
echo " * sample N : $sample_n "
echo " * BAM N : $bam_n "
echo " * INTERVALS: $bed "
echo " * DBSNP VCF: $dbsnp_vcf "
echo " * COSMIC VCF: $cosmic_vcf "
echo " * VCF original: $vcf_original "
echo " * VCF fixed: $vcf_fixed "

# single-threaded (multi-threading with -nct was much slower)

# --max_alt_allele_in_normal_fraction - threshold for maximum alternate allele fraction in normal [0.03]
# --max_alt_alleles_in_normal_count - threshold for maximum alternate allele counts in normal [1]
# --max_alt_alleles_in_normal_qscore_sum - threshold for maximum alternate allele quality score sum in normal [20]

mutect_cmd="
$gatk_cmd -T MuTect2 -dt NONE $gatk_log_level_arg \
--standard_min_confidence_threshold_for_calling 30 \
--max_alt_alleles_in_normal_count 10 \
--max_alt_allele_in_normal_fraction 0.05 \
--max_alt_alleles_in_normal_qscore_sum 40 \
--reference_sequence $ref_fasta \
$mt_var_arg \
--intervals $bed \
--interval_padding 10 \
--input_file:tumor $bam_t \
--input_file:normal $bam_n \
--out $vcf_original
"
echo -e "\n CMD: $mutect_cmd \n"
$mutect_cmd


#########################


# check that output generated

# check if VCF file is present
if [ ! -s "$vcf_original" ] ; then
	echo -e "\n $script_name ERROR: VCF $vcf_original NOT GENERATED \n" >&2
	exit 1
fi

# check if VCF index is present (should be present if VCF is complete)
if [ ! -s "$idx_original" ] ; then
	echo -e "\n $script_name ERROR: VCF IDX $idx_original NOT GENERATED \n" >&2
	# delete VCF since something went wrong and it might be corrupted
	rm -fv "$vcf_original"
	exit 1
fi


#########################


# adjust VCF for ANNOVAR compatibility (http://annovar.openbioinformatics.org/en/latest/articles/VCF/)

module unload samtools
module load samtools/1.3

# 1) keep header and only passing variants
# 2) split multi-allelic variants calls into separate lines (uses VCF 4.2 specification)
# 3) perform indel left-normalization (start position shifted to the left until it is no longer possible to do so)

fix_vcf_cmd="
cat $vcf_original \
| grep -E '^#|PASS' \
| bcftools norm --multiallelics -both --output-type v - \
| bcftools norm --fasta-ref $ref_fasta --output-type v - \
> $vcf_fixed
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
