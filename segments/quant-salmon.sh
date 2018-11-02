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

salmon_proj_logs_dir="${proj_dir}/logs-Salmon"
mkdir -p "$salmon_proj_logs_dir"
salmon_logs_dir="${salmon_proj_logs_dir}/${sample}"
salmon_quant_sf="${salmon_logs_dir}/quant.sf"
salmon_quant_genes_sf="${salmon_logs_dir}/quant.genes.sf"
salmon_quant_log="${salmon_logs_dir}/logs/salmon_quant.log"
salmon_quant_json="${salmon_logs_dir}/lib_format_counts.json"

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
conda activate salmon-0.11.3

#  -l [ --libType ] arg       Format string describing the library type
#  -i [ --index ] arg         salmon index
#  -r [ --unmatedReads ] arg  List of files containing unmated reads of (e.g. single-end reads)
#  -1 [ --mates1 ] arg        File containing the #1 mates
#  -2 [ --mates2 ] arg        File containing the #2 mates
#  --seqBias                  Perform sequence-specific bias correction.
#  --gcBias                   [beta for single-end reads] Perform fragment GC bias correction
#  -g [ --geneMap ] arg       File containing a mapping of transcripts to genes.  If this file is
#                             provided salmon will output both quant.sf and quant.genes.sf files,
#                             where the latter contains aggregated gene-level abundance estimates.
#  -q [ --quiet ]             Be quiet while doing quantification (no log file generated either)

# adjust for single/paired end
if [ -n "$fastq_R2" ] ;then
	fastq_param="--mates1 $fastq_R1 --mates2 $fastq_R2"
else
	fastq_param="--unmatedReads $fastq_R1"
fi

echo
echo " * conda version: $(conda --version | head -1) "
echo " * Salmon: $(readlink -f $(which salmon)) "
echo " * Salmon version: $(salmon --version | head -1) "
echo " * Salmon index: $salmon_index "
echo " * GTF: $gtf "
echo " * FASTQ R1: $fastq_R1 "
echo " * FASTQ R2: $fastq_R2 "
echo " * out dir: $salmon_logs_dir "
echo " * quant.genes.sf: $salmon_quant_genes_sf "
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

if [ ! -s "$salmon_quant_genes_sf" ] ; then
	echo -e "\n $script_name ERROR: RESULTS FILE $salmon_quant_genes_sf NOT GENERATED \n" >&2
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


# clean up

# clean up the full output table
# Name	Length	EffectiveLength	TPM	NumReads
echo -e "#GENE\t${sample}" > "$salmon_counts_txt"
cat "$salmon_quant_genes_sf" | grep -v "EffectiveLength" | cut -f 1,5 | LC_ALL=C sort -k1,1 >> "$salmon_counts_txt"
echo -e "#GENE\t${sample}" > "$salmon_tpms_txt"
cat "$salmon_quant_genes_sf" | grep -v "EffectiveLength" | cut -f 1,4 | LC_ALL=C sort -k1,1 >> "$salmon_tpms_txt"

# rename and gzip the quant.sf quantification file
mv -v "$salmon_quant_sf" "${salmon_quant_dir}/${sample}.quant.sf"
gzip "${salmon_quant_dir}/${sample}.quant.sf"


#########################


# generate summary

# jq binary (command-line JSON processor)
jq="/ifs/home/id460/software/jq/jq-1.5-linux64"

# extract stats
# lib_type=$(cat "$salmon_quant_log" | grep -m 1 "most likely library type" | sed 's/.*library type as //')
lib_type=$(cat "$salmon_quant_json" | "$jq" -r ".expected_format")
echo "library type: $lib_type"
# num_compatible_frags=$(cat "$salmon_quant_json" | "$jq" -r ".num_compatible_fragments")
# echo "compatible fragments: $num_compatible_frags"
num_assigned_frags=$(cat "$salmon_quant_json" | "$jq" -r ".num_assigned_fragments")
echo "assigned fragments: $num_assigned_frags"
map_rate=$(cat "$salmon_quant_log" | grep -m 1 "Mapping rate" | sed 's/.*rate = //')
echo "mapping rate: $map_rate"
num_genes=$(cat "$salmon_counts_txt" | grep -v "#GENE" | awk -F $'\t' '$2 > 0' | wc -l)
echo "detected genes: $num_genes"

# header for summary file
echo "#SAMPLE,LIBRARY TYPE, NUM ASSIGNED FRAGMENTS, MAPPING RATE, NUM GENES" > "$summary_csv"

# summarize log file
echo "${sample},${lib_type},${num_assigned_frags},${map_rate},${num_genes}" >> "$summary_csv"

sleep 30

# combine all sample summaries
cat ${summary_dir}/*.${segment_name}.csv | LC_ALL=C sort -t ',' -k1,1 | uniq > "${proj_dir}/summary.${segment_name}.csv"


#########################


# generate counts matrix for all samples
# this may be possible with "salmon quantmerge" in the future (does not currently support gene-level output)
# "tximport is, and has been, the recommended way to aggregate transcript-level abundances to the gene-level"

# merge counts if many samples are not still processing ("aux_info" directory is deleted at the end of this segment)
num_active_samples=$(find "$salmon_proj_logs_dir" -type d -name "aux_info" | wc -l)
if [ "$num_active_samples" -lt 5 ] ; then

	# load relevant modules
	module load r/3.3.0

	echo
	echo " * R: $(readlink -f $(which R)) "
	echo " * R version: $(R --version | head -1) "
	echo " * Rscript: $(readlink -f $(which Rscript)) "
	echo " * Rscript version: $(Rscript --version 2>&1) "
	echo

	Rscript --vanilla ${code_dir}/scripts/test-package.R optparse
	Rscript --vanilla ${code_dir}/scripts/test-package.R mnormt
	Rscript --vanilla ${code_dir}/scripts/test-package.R limma

	merged_counts_base="${proj_dir}/quant.salmon"

	# launch the analysis R script
	bash_cmd="Rscript --vanilla ${code_dir}/scripts/quant-merge-salmon.R $gtf $salmon_quant_dir $merged_counts_base"
	echo "CMD: $bash_cmd"
	($bash_cmd)

else

	echo
	echo "skip merging counts (quant-merge-salmon.R) since many samples are still processing"
	echo

fi


#########################


# clean up

gzip "${salmon_logs_dir}/quant.genes.sf"
rm -rfv "${salmon_logs_dir}/aux_info"
rm -rfv "${salmon_logs_dir}/libParams"


#########################



# end
