#!/bin/bash


# GATK coverage stats


# script filename
script_path="${BASH_SOURCE[0]}"
script_name=$(basename "$script_path")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 3 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir sample_name BAM \n" >&2
	if [ $# -gt 0 ] ; then echo -e "\n ARGS: $* \n" >&2 ; fi
	exit 1
fi

# arguments
proj_dir=$(readlink -f "$1")
sample=$2
bam=$3


#########################


# settings and files

summary_dir="${proj_dir}/summary"
mkdir -p "$summary_dir"
summary_csv="${summary_dir}/${sample}.${segment_name}.csv"

cov_dir="${proj_dir}/QC-coverage"
mkdir -p "$cov_dir"
out_prefix="${cov_dir}/${sample}"
gatk_sample_summary="${out_prefix}.sample_summary"

# unload all loaded modulefiles
module purge
module add default-environment


#########################


# exit if output exits already

if [ -s "$gatk_sample_summary" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 0
fi


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

ref_fasta=$(bash "${code_dir}/scripts/get-set-setting.sh" "${proj_dir}/settings.txt" REF-FASTA);

if [ ! -s "$ref_fasta" ] ; then
	echo -e "\n $script_name ERROR: FASTA $ref_fasta DOES NOT EXIST \n" >&2
	exit 1
fi

ref_dict=$(bash "${code_dir}/scripts/get-set-setting.sh" "${proj_dir}/settings.txt" REF-DICT);

if [ ! -s "$ref_dict" ] ; then
	echo -e "\n $script_name ERROR: DICT $ref_dict DOES NOT EXIST \n" >&2
	exit 1
fi

found_bed=$(find "$proj_dir" -maxdepth 1 -type f -iname "*.bed" | grep -v "probes" | sort | head -1)
bed=$(bash "${code_dir}/scripts/get-set-setting.sh" "${proj_dir}/settings.txt" EXP-TARGETS-BED "$found_bed");

if [ ! -s "$bed" ] ; then
	echo -e "\n $script_name ERROR: BED $bed DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# GATK settings

# command
# all other GATK segments work with 16G (hg19/mm10 WGS/WES)
# this segment failed for canFam3 WES (1.1M targets) with error "adjust the maximum heap size provided to Java"
gatk_jar="/gpfs/data/igorlab/software/GenomeAnalysisTK/GenomeAnalysisTK-3.8-1/GenomeAnalysisTK.jar"
gatk_cmd="java -Xms16G -Xmx16G -jar ${gatk_jar}"

if [ ! -s "$gatk_jar" ] ; then
	echo -e "\n $script_name ERROR: GATK $gatk_jar DOES NOT EXIST \n" >&2
	exit 1
fi

# error log (blank for troubleshooting)
gatk_log_level_arg="--logging_level ERROR"


#########################


# on-target coverage

echo
echo " * GATK: $(readlink -f $gatk_jar) "
echo " * GATK version: $($gatk_cmd --version) "
echo " * BAM: $bam "
echo " * BED: $bed "
echo " * output prefix: $out_prefix "
echo " * sample_summary: $gatk_sample_summary "
echo

# using '-nt' with this combination of arguments causes an error

gatk_doc_cmd="
$gatk_cmd -T DepthOfCoverage -dt NONE $gatk_log_level_arg \
-rf BadCigar \
--reference_sequence $ref_fasta \
--intervals $bed \
--omitDepthOutputAtEachBase \
-ct 10 -ct 50 -ct 100 -ct 500 -mbq 20 -mmq 20 --nBins 999 --start 1 --stop 1000 \
--input_file $bam \
--outputFormat csv \
--out $out_prefix
"
echo "CMD: $gatk_doc_cmd"
$gatk_doc_cmd


#########################


# check that output generated

if [ ! -s "$gatk_sample_summary" ] ; then
	echo -e "\n $script_name ERROR: sample_summary $gatk_sample_summary NOT GENERATED \n" >&2
	exit 1
fi


#########################


# generate summary

# summarize log file
cat "$gatk_sample_summary" \
| head -2 \
| cut -d ',' -f 1,3,5,7-99 \
| sed 's/sample_id,mean,granular_median/#SAMPLE,MEAN COVERAGE,MEDIAN COVERAGE/' \
> "$summary_csv"

sleep 5

# combine all sample summaries
cat ${summary_dir}/*.${segment_name}.csv | LC_ALL=C sort -t ',' -k1,1 | uniq > "${proj_dir}/summary.${segment_name}.csv"


#########################


# delete files that are not needed
rm -fv "${out_prefix}.sample_statistics"
rm -fv "${out_prefix}.sample_interval_statistics"


#########################



# end
