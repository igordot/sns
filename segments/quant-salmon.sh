#!/bin/bash


# run Salmon (quasi-mapping-based mode)


# script filename
script_path="${BASH_SOURCE[0]}"
script_name=$(basename "$script_path")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ $# -lt 4 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir sample_name num_threads FASTQ_R1 [FASTQ_R2] \n" >&2
	exit 1
fi

# arguments
proj_dir=$1
sample=$2
threads=$3
fastq_R1=$4
fastq_R2=$5


#########################


# settings and files

summary_dir="${proj_dir}/summary"
mkdir -p "$summary_dir"
summary_csv="${summary_dir}/${sample}.${segment_name}.csv"

salmon_quant_dir="${proj_dir}/quant-Salmon"
mkdir -p "$salmon_quant_dir"
salmon_counts_txt="${salmon_quant_dir}/${sample}.reads.txt"
salmon_tpms_txt="${salmon_quant_dir}/${sample}.tpms.txt"

salmon_logs_dir="${proj_dir}/logs-Salmon"
mkdir -p "$salmon_logs_dir"
salmon_logs_dir="${salmon_logs_dir}/${sample}"
salmon_quant_sf="${salmon_logs_dir}/quant.genes.sf"
salmon_quant_log="${salmon_logs_dir}/logs/salmon_quant.log"

# unload all loaded modulefiles
module purge
module load local


#########################


# check for output

if [ -s "$salmon_counts_txt" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 1
fi


#########################


# check that inputs exist

if [ ! -d "$proj_dir" ] ; then
	echo -e "\n $script_name ERROR: DIR $proj_dir DOES NOT EXIST \n" >&2
	exit 1
fi

if [ ! -s "$fastq_R1" ] ; then
	echo -e "\n $script_name ERROR: FASTQ $fastq_R1 DOES NOT EXIST \n" >&2
	exit 1
fi

code_dir=$(dirname $(dirname "$script_path"))

gtf=$(bash "${code_dir}/scripts/get-set-setting.sh" "${proj_dir}/settings.txt" REF-GTF);

if [ ! -s "$gtf" ] || [ ! "$gtf" ] ; then
	echo -e "\n $script_name ERROR: GTF $gtf DOES NOT EXIST \n" >&2
	exit 1
fi

salmon_index=$(bash "${code_dir}/scripts/get-set-setting.sh" "${proj_dir}/settings.txt" REF-SALMON);

if [ ! -s "${salmon_index}/hash.bin" ] ; then
	echo -e "\n $script_name ERROR: SALMON INDEX $salmon_index DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# Salmon

# prepare and activate Conda environment
module purge
module load local
unset PYTHONPATH
source /ifs/home/id460/anaconda3/etc/profile.d/conda.sh
conda activate salmon-0.10.2

# adjust for single/paired end
if [ -n "$fastq_R2" ] ;then
	fastq_param="--mates1 $fastq_R1 --mates2 $fastq_R2"
else
	fastq_param="--unmatedReads $fastq_R1"
fi

echo
echo " * conda: $(readlink -f $(which conda)) "
echo " * conda version: $(conda --version | head -1) "
echo " * Salmon: $(readlink -f $(which salmon)) "
echo " * Salmon version: $(salmon --version | head -1) "
echo " * Salmon index: $salmon_index "
echo " * GTF: $gtf "
echo " * FASTQ R1: $fastq_R1 "
echo " * FASTQ R2: $fastq_R2 "
echo " * out dir: $salmon_logs_dir "
echo " * quant.genes.sf: $salmon_quant_sf "
echo " * reads: $salmon_counts_txt "
echo " * TPMs: $salmon_tpms_txt "
echo

salmon_cmd="
salmon --no-version-check quant \
--threads $threads \
--index $salmon_index \
--geneMap $gtf \
--libType A \
--seqBias --gcBias \
$fastq_param \
--output $salmon_logs_dir
"
echo -e "\n CMD: $salmon_cmd \n"
$salmon_cmd

# deactivate Conda environment
conda deactivate

sleep 30


#########################


# check that output generated

if [ ! -s "$salmon_quant_sf" ] ; then
	echo -e "\n $script_name ERROR: RESULTS FILE $salmon_quant_sf NOT GENERATED \n" >&2
	# delete supplementary files since something went wrong and they might be corrupted
	rm -rfv "${salmon_logs_dir}/aux_info"
	rm -rfv "${salmon_logs_dir}/libParams"
	exit 1
fi

if [ ! -s "$salmon_quant_log" ] ; then
	echo -e "\n $script_name ERROR: LOG FILE $salmon_quant_log NOT GENERATED \n" >&2
	# delete supplementary files since something went wrong and they might be corrupted
	rm -rfv "${salmon_logs_dir}/aux_info"
	rm -rfv "${salmon_logs_dir}/libParams"
	exit 1
fi


#########################


# clean up the full output table

# Name	Length	EffectiveLength	TPM	NumReads

echo -e "#GENE\t${sample}" > "$salmon_counts_txt"
cat "$salmon_quant_sf" | grep -v "EffectiveLength" | cut -f 1,5 | LC_ALL=C sort -k1,1 >> "$salmon_counts_txt"

echo -e "#GENE\t${sample}" > "$salmon_tpms_txt"
cat "$salmon_quant_sf" | grep -v "EffectiveLength" | cut -f 1,4 | LC_ALL=C sort -k1,1 >> "$salmon_tpms_txt"


#########################


# generate summary

# extract stats (also consider parsing lib_format_counts.json)
library_type=$(cat "$salmon_quant_log" | grep -m 1 "most likely library type" | sed 's/.*library type as //')
echo "library type: $library_type"
map_rate=$(cat "$salmon_quant_log" | grep -m 1 "Mapping rate" | sed 's/.*rate = //')
echo "mapping rate: $map_rate"

# header for summary file
echo "#SAMPLE,LIBRARY TYPE,MAPPING RATE" > "$summary_csv"

# summarize log file
echo "${sample},${library_type},${map_rate}" >> "$summary_csv"

sleep 30

# combine all sample summaries
cat ${summary_dir}/*.${segment_name}.csv | LC_ALL=C sort -t ',' -k1,1 | uniq > "${proj_dir}/summary.${segment_name}.csv"


#########################


# generate counts matrix for all samples
# this may be possible with "salmon quantmerge" in the future (does not currently support gene-level output)
# add sample to temp merged filename to avoid conflicts during join process

merged_reads_txt="${proj_dir}/quant.salmon.reads.txt"
merged_tpms_txt="${proj_dir}/quant.salmon.tpms.txt"

reads_file_list=$(find "$salmon_quant_dir" -type f -name "*.reads.txt" | LC_ALL=C sort | tr '\n' ' ')
join_cmd="bash ${code_dir}/scripts/join-many.sh $'\t' 0 $reads_file_list > ${merged_reads_txt}.${sample}.tmp"
(eval "$join_cmd")
mv -fv "${merged_reads_txt}.${sample}.tmp" "$merged_reads_txt"

tpms_file_list=$(find "$salmon_quant_dir" -type f -name "*.tpms.txt" | LC_ALL=C sort | tr '\n' ' ')
join_cmd="bash ${code_dir}/scripts/join-many.sh $'\t' 0 $tpms_file_list > ${merged_tpms_txt}.${sample}.tmp"
(eval "$join_cmd")
mv -fv "${merged_tpms_txt}.${sample}.tmp" "$merged_tpms_txt"


#########################


# clean up

gzip "${salmon_logs_dir}/quant.sf"
gzip "${salmon_logs_dir}/quant.genes.sf"
rm -rfv "${salmon_logs_dir}/aux_info"
rm -rfv "${salmon_logs_dir}/libParams"


#########################



# end
