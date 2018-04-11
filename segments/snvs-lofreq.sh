#!/bin/bash


# call variants with LoFreq


# script filename
script_path="${BASH_SOURCE[0]}"
script_name=$(basename "$script_path")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 4 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir sample_name threads BAM \n" >&2
	if [ $# -gt 0 ] ; then echo -e "\n ARGS: $* \n" >&2 ; fi
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

code_dir=$(dirname $(dirname "$script_path"))

ref_fasta=$(bash "${code_dir}/scripts/get-set-setting.sh" "${proj_dir}/settings.txt" REF-FASTA)

if [ ! -s "$ref_fasta" ] ; then
	echo -e "\n $script_name ERROR: FASTA $ref_fasta DOES NOT EXIST \n" >&2
	exit 1
fi

chrom_sizes=$(bash "${code_dir}/scripts/get-set-setting.sh" "${proj_dir}/settings.txt" REF-CHROMSIZES)

if [ ! -s "$chrom_sizes" ] ; then
	echo -e "\n $script_name ERROR: CHROM SIZES $chrom_sizes DOES NOT EXIST \n" >&2
	exit 1
fi

found_bed=$(find "$proj_dir" -maxdepth 1 -type f -iname "*.bed" | grep -v "probes" | sort | head -1)
bed=$(bash "${code_dir}/scripts/get-set-setting.sh" "${proj_dir}/settings.txt" EXP-TARGETS-BED "$found_bed")

if [ ! -s "$bed" ] ; then
	echo -e "\n $script_name ERROR: BED $bed DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# settings and files

vcf_dir="${proj_dir}/VCF-LoFreq"
mkdir -p "$vcf_dir"
vcf_original="${vcf_dir}/${sample}.original.vcf"
vcf_add_gt="${vcf_dir}/${sample}.gt.vcf"
vcf_fixed="${vcf_dir}/${sample}.vcf"

# annotation command (next segment)
annot_cmd="bash ${code_dir}/segments/annot-annovar.sh $proj_dir $sample $vcf_fixed"

# unload all loaded modulefiles
module purge
module load local


#########################


# check for output

# skip to annotation if final output exists already
if [ -s "$vcf_fixed" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	echo -e "\n CMD: $annot_cmd \n"
	($annot_cmd)
	exit 1
fi

# delete original VCF (likely incomplete since the fixed VCF was not generated)
if [ -s "$vcf_original" ] ; then
	echo -e "\n $script_name WARNING: POTENTIALLY CORRUPT VCF $vcf_original EXISTS \n" >&2
	rm -fv "$vcf_original"
fi


#########################


# create padded bed if needed

# two steps in case orginal BED file does not end in ".bed"
bed_padded="${bed}.pad10"
bed_padded=${bed_padded/%.bed.pad10/.pad10.bed}

bed_pad_cmd="
cat $bed \
| LC_ALL=C sort -k1,1 -k2,2n \
| bedtools slop -g $chrom_sizes -b 10 \
| bedtools merge -d 5 \
> $bed_padded
"

if [ ! -s "$bed_padded" ] ; then
	module load bedtools/2.26.0
	echo -e "\n CMD: $bed_pad_cmd \n"
	eval "$bed_pad_cmd"
	sleep 30
fi

if [ ! -s "$bed_padded" ] ; then
	echo -e "\n $script_name ERROR: BED $bed_padded DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# LoFreq

lofreq_bin="/ifs/home/id460/bin/lofreq"
lofreq_gt_py="/ifs/home/id460/bin/lofreq2_add_fake_gt.py"

echo
echo " * LoFreq: $(readlink -f $(which $lofreq_bin)) "
echo " * LoFreq version: $($lofreq_bin version 2>&1 | head -1) "
echo " * BAM: $bam "
echo " * BED: $bed_padded "
echo " * VCF original: $vcf_original "
echo " * VCF fixed: $vcf_fixed "
echo

lofreq_cmd="
$lofreq_bin call-parallel --call-indels \
--pp-threads $threads \
--ref $ref_fasta \
--bed $bed_padded \
--out $vcf_original \
$bam
"
echo -e "\n CMD: $lofreq_cmd \n"
$lofreq_cmd


#########################


# check that output generated

if [ ! -s "$vcf_original" ] ; then
	echo -e "\n $script_name ERROR: VCF $vcf_original NOT GENERATED \n" >&2
	exit 1
fi


#########################


# complement VCF with unknown genotype (adds FORMAT and SAMPLE fields for improved compatibility)
# https://github.com/CSB5/lofreq/blob/master/src/tools/scripts/lofreq2_add_fake_gt.py
# "you can just add fake columns"
# http://csb5.github.io/lofreq/2015/11/23/where-are-the-format-and-sample-fields/

module load python/2.7.3

echo
echo " * lofreq2_add_fake_gt.py: $(readlink -f $(which $lofreq_gt_py)) "
echo " * VCF original: $vcf_original "
echo " * VCF with GT: $vcf_add_gt "
echo

lofreq_gt_cmd="
python $lofreq_gt_py \
--vcf-in $vcf_original \
--samples $sample \
--vcf-out $vcf_add_gt \
"
echo -e "\n CMD: $lofreq_gt_cmd \n"
$lofreq_gt_cmd


#########################


# check that output generated

if [ ! -s "$vcf_add_gt" ] ; then
	echo -e "\n $script_name ERROR: VCF $vcf_add_gt NOT GENERATED \n" >&2
	exit 1
fi


#########################


# adjust the vcf for annovar compatibility (http://www.openbioinformatics.org/annovar/annovar_vcf.html)

module load samtools/1.3

echo
echo " * samtools: $(readlink -f $(which samtools)) "
echo " * samtools version: $(samtools --version | head -1) "
echo

# create indexed VCF file
# warnings: "contig 'chr1' is not defined in the header. (Quick workaround: index the file with tabix.)"

# bgzip VCF file so it can be indexed
bgzip_cmd="bgzip -c $vcf_add_gt > ${vcf_original}.bgz"
echo -e "\n CMD: $bgzip_cmd \n"
eval "$bgzip_cmd"

# index VCF file
bcf_index_cmd="bcftools index ${vcf_original}.bgz"
echo -e "\n CMD: $bcf_index_cmd \n"
eval "$bcf_index_cmd"

# 1) split multi-allelic variants calls into separate lines (uses VCF 4.2 specification)
# 2) perform indel left-normalization (start position shifted to the left until it is no longer possible to do so)
# 3) depth filter

fix_vcf_cmd="
bcftools norm --multiallelics -both --output-type v ${vcf_original}.bgz \
| bcftools norm --fasta-ref $ref_fasta --output-type v - \
| bcftools view --exclude 'DP<5' --output-type v > $vcf_fixed
"
echo -e "\n CMD: $fix_vcf_cmd \n"
eval "$fix_vcf_cmd"


#########################


# check that output generated

if [ ! -s "$vcf_fixed" ] ; then
	echo -e "\n $script_name ERROR: VCF $vcf_fixed NOT GENERATED \n" >&2
	# delete all VCFs since something went wrong and they might be corrupted
	rm -fv "$vcf_original"
	rm -fv "${vcf_original}.bgz"
	rm -fv "${vcf_original}.bgz.csi"
	rm -fv "$vcf_fixed"
	exit 1
fi


#########################


# clean up

rm -fv "$vcf_add_gt"
rm -fv "${vcf_original}.bgz"
rm -fv "${vcf_original}.bgz.csi"


#########################


# annotate

echo -e "\n CMD: $annot_cmd \n"
($annot_cmd)


#########################



# end
