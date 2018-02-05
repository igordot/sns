#!/bin/bash


# call somatic variants with Strelka


# script filename
script_path="${BASH_SOURCE[0]}"
script_name=$(basename "$script_path")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 6 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir tumor_sample_name tumor_bam normal_sample_name normal_bam threads \n" >&2
	exit 1
fi

# arguments
proj_dir=$1
sample_t=$2
bam_t=$3
sample_n=$4
bam_n=$5
threads=$6


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

code_dir=$(dirname $(dirname "$script_path"))

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


#########################


# settings and files

sample="${sample_t}:${sample_n}"

# Manta output directory (best-practice recommendation is to run the Manta SV and indel caller first)
manta_logs_dir="${proj_dir}/logs-${segment_name}/${sample_t}-${sample_n}-manta"

manta_run_py="${manta_logs_dir}/runWorkflow.py"
manta_indels_vcf="${manta_logs_dir}/results/variants/candidateSmallIndels.vcf.gz"

# Strelka output directory (must not exist before running Strelka)
strelka_logs_dir="${proj_dir}/logs-${segment_name}/${sample_t}-${sample_n}-strelka"

strelka_run_py="${strelka_logs_dir}/runWorkflow.py"

vcf_dir="${proj_dir}/VCF-Strelka"
mkdir -p "$vcf_dir"
vcf_snvs_original="${vcf_dir}/${sample_t}-${sample_n}.snvs.original.vcf.gz"
vcf_indels_original="${vcf_dir}/${sample_t}-${sample_n}.indels.original.vcf.gz"
vcf_snvs_fixed="${vcf_dir}/${sample_t}-${sample_n}.snvs.vcf"
vcf_indels_fixed="${vcf_dir}/${sample_t}-${sample_n}.indels.vcf"
vcf_combined="${vcf_dir}/${sample_t}-${sample_n}.all.vcf"

# annotation command (next segment)
annot_cmd="bash ${code_dir}/segments/annot-annovar.sh $proj_dir $sample $vcf_combined"

# unload all loaded modulefiles
module purge
module load local


#########################


# check for output

# skip to annotation if final output exists already
if [ -s "$vcf_combined" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	echo -e "\n CMD: $annot_cmd \n"
	($annot_cmd)
	exit 1
fi

# delete candidateSmallIndels.vcf.gz (likely incomplete since the final VCF was not generated)
if [ -s "$manta_indels_vcf" ] ; then
	echo -e "\n $script_name WARNING: POTENTIALLY CORRUPT VCF $manta_indels_vcf EXISTS \n" >&2
	rm -fv "$manta_indels_vcf"
fi

# delete original SNVs VCF (likely incomplete since the final VCF was not generated)
if [ -s "$vcf_snvs_original" ] ; then
	echo -e "\n $script_name WARNING: POTENTIALLY CORRUPT VCF $vcf_snvs_original EXISTS \n" >&2
	rm -fv "$vcf_snvs_original"
fi

# delete original indels VCF (likely incomplete since the final VCF was not generated)
if [ -s "$vcf_indels_original" ] ; then
	echo -e "\n $script_name WARNING: POTENTIALLY CORRUPT VCF $vcf_indels_original EXISTS \n" >&2
	rm -fv "$vcf_indels_original"
fi


#########################


# configure Manta

# Manta/Strelka use python-based config scripts
module load python/2.7.3

manta_dir="/ifs/home/id460/software/manta/manta-1.3.0"
manta_config_py="${manta_dir}/bin/configManta.py"

echo
echo " * Manta config: $(readlink -f $manta_config_py) "
echo " * Manta version: $($manta_config_py --version) "
echo " * sample T : $sample_t "
echo " * BAM T : $bam_t "
echo " * sample N : $sample_n "
echo " * BAM N : $bam_n "
echo " * output dir : $manta_logs_dir "
echo " * indels VCF : $manta_indels_vcf "
echo

manta_config_cmd="
$manta_config_py \
--exome \
--referenceFasta $ref_fasta \
--tumorBam $bam_t \
--normalBam $bam_n \
--runDir $manta_logs_dir
"
echo -e "\n CMD: $manta_config_cmd \n"
$manta_config_cmd


#########################


# check that output generated

# check if runWorkflow.py is present
if [ ! -s "$manta_run_py" ] ; then
	echo -e "\n $script_name ERROR: runWorkflow.py $manta_run_py NOT GENERATED \n" >&2
	# delete Manta output (keep top level for logs)
	rm -rfv "${manta_logs_dir}/results"
	rm -rfv "${manta_logs_dir}/workspace"
	exit 1
fi


#########################


# run Manta

manta_run_cmd="
$manta_run_py \
--quiet \
--mode local \
--jobs $threads \
--memGb 64 \
"
echo -e "\n CMD: $manta_run_cmd \n"
$manta_run_cmd

sleep 30


#########################


# check that output generated

# check if candidateSmallIndels.vcf.gz is present
if [ ! -s "$manta_indels_vcf" ] ; then
	echo -e "\n $script_name ERROR: VCF $manta_indels_vcf NOT GENERATED \n" >&2
	# delete Manta output (keep top level for logs)
	rm -rfv "${manta_logs_dir}/results"
	rm -rfv "${manta_logs_dir}/workspace"
	exit 1
fi


#########################


# configure Strelka

strelka_dir="/ifs/home/id460/software/strelka/strelka-2.8.4"
strelka_config_py="${strelka_dir}/bin/configureStrelkaSomaticWorkflow.py"

echo
echo " * Strelka config: $(readlink -f $strelka_config_py) "
echo " * Strelka version: $($strelka_config_py --version) "
echo " * sample T : $sample_t "
echo " * BAM T : $bam_t "
echo " * sample N : $sample_n "
echo " * BAM N : $bam_n "
echo " * output dir : $strelka_logs_dir "
echo " * VCF original: $vcf_snvs_original "
echo " * VCF original: $vcf_indels_original "
echo " * VCF fixed: $vcf_snvs_fixed "
echo " * VCF fixed: $vcf_indels_fixed "
echo

strelka_config_cmd="
$strelka_config_py \
--exome \
--outputCallableRegions \
--referenceFasta $ref_fasta \
--tumorBam $bam_t \
--normalBam $bam_n \
--indelCandidates $manta_indels_vcf \
--runDir $strelka_logs_dir
"
echo -e "\n CMD: $strelka_config_cmd \n"
$strelka_config_cmd


#########################


# check that output generated

# check if runWorkflow.py is present
if [ ! -s "$strelka_run_py" ] ; then
	echo -e "\n $script_name ERROR: runWorkflow.py $strelka_run_py NOT GENERATED \n" >&2
	# delete Manta and Strelka output (keep top level for logs)
	rm -rfv "${manta_logs_dir}"
	rm -rfv "${strelka_logs_dir}/results"
	rm -rfv "${strelka_logs_dir}/workspace"
	exit 1
fi


#########################


# run Strelka

strelka_run_cmd="
$strelka_run_py \
--quiet \
--mode local \
--jobs $threads \
--memGb 64 \
"
echo -e "\n CMD: $strelka_run_cmd \n"
$strelka_run_cmd

sleep 30

# move and rename the default VCFs
mv -v "${strelka_logs_dir}/results/variants/somatic.snvs.vcf.gz" "$vcf_snvs_original"
mv -v "${strelka_logs_dir}/results/variants/somatic.indels.vcf.gz" "$vcf_indels_original"

sleep 30


#########################


# check that output generated

# check if SNVs VCF file is present
if [ ! -s "$vcf_snvs_original" ] ; then
	echo -e "\n $script_name ERROR: VCF $vcf_snvs_original NOT GENERATED \n" >&2
	# delete Manta and Strelka output (keep top level for logs)
	rm -rfv "${manta_logs_dir}"
	rm -rfv "${strelka_logs_dir}/results"
	rm -rfv "${strelka_logs_dir}/workspace"
	exit 1
fi

# check if indels VCF file is present
if [ ! -s "$vcf_indels_original" ] ; then
	echo -e "\n $script_name ERROR: VCF $vcf_indels_original NOT GENERATED \n" >&2
	# delete Manta and Strelka output (keep top level for logs)
	rm -rfv "${manta_logs_dir}"
	rm -rfv "${strelka_logs_dir}/results"
	rm -rfv "${strelka_logs_dir}/workspace"
	exit 1
fi


#########################


# clean up

rm -rf "${manta_logs_dir}/workspace"
rm -rf "${strelka_logs_dir}/workspace"


#########################


# adjust VCF for ANNOVAR compatibility (http://annovar.openbioinformatics.org/en/latest/articles/VCF/)

module load samtools/1.3

# 1) keep header and only passing variants
# 2) split multi-allelic variants calls into separate lines (uses VCF 4.2 specification)
# 3) perform indel left-normalization (start position shifted to the left until it is no longer possible to do so)

fix_snvs_vcf_cmd="
gzip -cd $vcf_snvs_original \
| grep -E '^#|PASS' \
| bcftools norm --multiallelics -both --output-type v - \
| bcftools norm --fasta-ref $ref_fasta --output-type v - \
> $vcf_snvs_fixed
"
echo -e "\n CMD: $fix_snvs_vcf_cmd \n"
eval "$fix_snvs_vcf_cmd"

fix_indels_vcf_cmd="
gzip -cd $vcf_indels_original \
| grep -E '^#|PASS' \
| bcftools norm --multiallelics -both --output-type v - \
| bcftools norm --fasta-ref $ref_fasta --output-type v - \
> $vcf_indels_fixed
"
echo -e "\n CMD: $fix_indels_vcf_cmd \n"
eval "$fix_indels_vcf_cmd"

sleep 30


#########################


# check that output generated

if [ ! -s "$vcf_snvs_fixed" ] ; then
	echo -e "\n $script_name ERROR: VCF $vcf_snvs_fixed NOT GENERATED \n" >&2
	exit 1
fi

if [ ! -s "$vcf_indels_fixed" ] ; then
	echo -e "\n $script_name ERROR: VCF $vcf_indels_fixed NOT GENERATED \n" >&2
	exit 1
fi


#########################


# create a combined VCF (produces a valid VCF according to bcftools)

cat \
<(cat "$vcf_snvs_fixed" | grep '^##' | grep -E -v '##content|##cmdline|##prior|##bcftools' ) \
<(cat "$vcf_indels_fixed" | grep '^##INFO') \
<(cat "$vcf_indels_fixed" | grep '^##FORMAT') \
<(cat "$vcf_indels_fixed" | grep '^#CHROM') \
<(cat "$vcf_snvs_fixed" "$vcf_indels_fixed" | grep -v '^#' | LC_ALL=C sort -k1,1 -k2,2n) \
> "$vcf_combined"

sleep 30


#########################


# check that output generated

if [ ! -s "$vcf_combined" ] ; then
	echo -e "\n $script_name ERROR: VCF $vcf_combined NOT GENERATED \n" >&2
	exit 1
fi


#########################


# annotate

echo -e "\n CMD: $annot_cmd \n"
($annot_cmd)


#########################



# end
