#!/bin/bash


# GATK coverage stats


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

found_bed=$(find $proj_dir -maxdepth 1 -type f -name "*.bed" | head -1)
bed=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" EXP-BED $found_bed);

if [ ! -s "$bed" ] ; then
	echo -e "\n $script_name ERROR: BED $bed DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# settings and files

summary_dir="${proj_dir}/summary"
mkdir -p "$summary_dir"
summary_csv="${summary_dir}/${sample}.${segment_name}.csv"

cov_dir="${proj_dir}/QC-coverage"
mkdir -p "$cov_dir"
out_prefix="${cov_dir}/${sample}"
gatk_sample_summary="${out_prefix}.sample_summary"

# logs_dir="${proj_dir}/logs-${segment_name}"
# mkdir -p "$logs_dir"


#########################


# exit if output exits already

if [ -s "$gatk_sample_summary" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 1
fi


#########################


# GATK settings

module unload java
module load java/1.8

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


# on-target coverage

echo " * GATK: $(readlink -f $gatk_jar) "
echo " * GATK version: $($gatk_cmd --version) "
echo " * BAM IN: $bam "
echo " * OUT PREFIX: $out_prefix "
echo " * sample_summary: $gatk_sample_summary "

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

sleep 30

# combine all sample summaries
cat ${summary_dir}/*.${segment_name}.csv | LC_ALL=C sort -t ',' -k1,1 | uniq > "${proj_dir}/summary.${segment_name}.csv"


#########################


# delete files that are not needed
rm -fv "${out_prefix}.sample_statistics"
rm -fv "${out_prefix}.sample_interval_statistics"


#########################



# end
