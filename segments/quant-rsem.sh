#!/bin/bash


# run RSEM


# script filename
script_name=$(basename "${BASH_SOURCE[0]}")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ $# -lt 5 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name [project dir] [sample] [threads] [strand] [FASTQ R1] [FASTQ R2] \n" >&2
	exit 1
fi

# arguments
proj_dir=$1
sample=$2
threads=$3
strand=$4
fastq_R1=$5
fastq_R2=$6


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

code_dir=$(dirname "$(dirname "${BASH_SOURCE[0]}")")
rsem_genome=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" REF-RSEM);

if [ ! -s "${rsem_genome}.transcripts.fa" ] ; then
	echo -e "\n $script_name ERROR: GENOME ${rsem_genome}.transcripts.fa DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# settings and files

summary_dir="${proj_dir}/summary"
mkdir -p "$summary_dir"
summary_csv="${summary_dir}/${sample}.rsem.${strand}.csv"

rsem_quant_dir="${proj_dir}/quant-RSEM"
mkdir -p "$rsem_quant_dir"
rsem_expcounts_txt="${rsem_quant_dir}/${sample}.expcounts.${strand}.txt"
rsem_tpm_txt="${rsem_quant_dir}/${sample}.TPM.${strand}.txt"
rsem_fpkm_txt="${rsem_quant_dir}/${sample}.FPKM.${strand}.txt"

rsem_logs_dir="${proj_dir}/logs-RSEM"
mkdir -p "$rsem_logs_dir"
# RSEM_RAW_SUMMARY="${rsem_logs_dir}/${sample}.${strand}.stats.txt"
rsem_raw_genes_results="${rsem_logs_dir}/${sample}.${strand}.genes.results"
rsem_stat_cnt="${rsem_logs_dir}/${sample}.${strand}.stat/${sample}.${strand}.cnt"

# adjust for single/paired end
if [ -n "$fastq_R2" ] ;then
	fastq_flag="--paired-end $fastq_R1 $fastq_R2"
else
	fastq_flag="$fastq_R1"
fi

# adjust for strand
# fwd | transcript             | cufflinks "fr-secondstrand" | htseq "yes"     | picard "FIRST_READ"
# rev | rev comp of transcript | cufflinks "fr-firststrand"  | htseq "reverse" | picard "SECOND_READ"
if [ "$strand" == "unstr" ] ; then
	strand_flag=""
elif [ "$strand" == "fwd" ] ; then
	strand_flag="--forward-prob 1"
elif [ "$strand" == "rev" ] ; then
	strand_flag="--forward-prob 0"
else
	echo -e "\n $script_name ERROR: INCORRECT strand $strand SPECIFIED \n" >&2
	exit 1
fi


#########################


# exit if output exits already

if [ -s "$rsem_expcounts_txt" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 1
fi


#########################


# RSEM

module load bowtie2/2.2.6

rsem_bin_dir="/ifs/home/id460/software/RSEM-1.2.31/bin"

echo " * RSEM: ${rsem_bin_dir}/rsem-calculate-expression "
echo " * RSEM version: $(${rsem_bin_dir}/rsem-calculate-expression --version) "
echo " * bowtie2: $(readlink -f $(which bowtie2)) "
echo " * bowtie2 version: $(bowtie2 --version 2>&1 | head -1) "
echo " * RSEM REF: $rsem_genome "
echo " * FASTQ R1: $fastq_R1 "
echo " * FASTQ R2: $fastq_R2 "
echo " * COUNTS : $rsem_expcounts_txt "
echo " * TPM : $rsem_tpm_txt "
echo " * FPKM : $rsem_fpkm_txt "

CMD="
${rsem_bin_dir}/rsem-calculate-expression \
--quiet \
--no-bam-output \
--bowtie2 \
--num-threads $threads \
--estimate-rspd \
$strand_flag \
$fastq_flag \
$rsem_genome \
${rsem_logs_dir}/${sample}.${strand}
"
echo "CMD: $CMD"
eval "$CMD"

sleep 30


#########################


# check that output generated

if [ ! -s "$rsem_raw_genes_results" ] ; then
	echo -e "\n $script_name ERROR: RESULTS FILE $rsem_raw_genes_results NOT GENERATED \n" >&2
	exit 1
fi


#########################


# clean up the full output

# all results
# cat $RSEM_GENES_RESULTS | sed 's/gene_id/000_gene_id/' | cut -f 1,5,6,7 | LC_ALL=C sort -k1,1 > $RSEM_GENES_RESULTS_TXT

echo -e "#GENE\t${sample}" > "$rsem_expcounts_txt"
cat "$rsem_raw_genes_results" | grep -v "gene_id" | cut -f 1,5 | LC_ALL=C sort -k1,1 >> "$rsem_expcounts_txt"

echo -e "#GENE\t${sample}" > "$rsem_tpm_txt"
cat "$rsem_raw_genes_results" | grep -v "gene_id" | cut -f 1,6 | LC_ALL=C sort -k1,1 >> "$rsem_tpm_txt"

echo -e "#GENE\t${sample}" > "$rsem_fpkm_txt"
cat "$rsem_raw_genes_results" | grep -v "gene_id" | cut -f 1,7 | LC_ALL=C sort -k1,1 >> "$rsem_fpkm_txt"


#########################


# generate summary

# *.cnt file contains alignment statistics based purely on the alignment results obtained from aligners
# N0 N1 N2 N_tot   #  N0, number of unalignable reads; N1, number of alignable reads; N2, number of filtered reads due to too many alignments; N_tot = N0 + N1 + N2
# nUnique nMulti nUncertain   # nUnique, number of reads aligned uniquely to a gene;
#                             # nMulti, number of reads aligned to multiple genes; nUnique + nMulti = N1;
#                             # nUncertain, number of reads aligned to multiple locations in the given reference sequences, which include isoform-level multi-mapping reads
# nHits read_type             # nHits, number of total alignments.
#                             # read_type: 0, single-end read, no quality score; 1, single-end read, with quality score; 2, paired-end read, no quality score; 3, paired-end read, with quality score

# extract read counts
reads_input=$(     cat "$rsem_stat_cnt" | head -1 | cut -d " " -f 4)
reads_map=$(       cat "$rsem_stat_cnt" | head -1 | cut -d " " -f 2)
reads_map_unique=$(cat "$rsem_stat_cnt" | head -2 | tail -1 | cut -d " " -f 1)
reads_map_multi=$( cat "$rsem_stat_cnt" | head -2 | tail -1 | cut -d " " -f 2)

# calculate ratios
map_pct=$(       echo "(${reads_map}/${reads_input})*100"        | bc -l | cut -c 1-4)
map_unique_pct=$(echo "(${reads_map_unique}/${reads_input})*100" | bc -l | cut -c 1-4)
map_multi_pct=$( echo "(${reads_map_multi}/${reads_input})*100"  | bc -l | cut -c 1-4)

map_pct="${map_pct}%"
map_unique_pct="${map_unique_pct}%"
map_multi_pct="${map_multi_pct}%"

# header for summary file
echo "#SAMPLE,INPUT READS,MAPPED READS,UNIQUELY MAPPED,MULTI-MAPPED,\
MAPPED READS,UNIQUELY MAPPED,MULTI-MAPPED" > "$summary_csv"

# summarize log file
echo "${sample},${reads_input},${reads_map},${reads_map_unique},${reads_map_multi},\
${map_pct},${map_unique_pct},${map_multi_pct}" >> "$summary_csv"

sleep 30

# combine all sample summaries
cat ${summary_dir}/*.rsem.${strand}.csv | LC_ALL=C sort -t ',' -k1,1 | uniq > "${proj_dir}/summary.rsem.${strand}.csv"


#########################


# generate counts matrix for all samples
# merged counts filename with sample included to avoid conflicts during join process

merged_counts_file="${proj_dir}/quant.rsem.expcounts.${strand}.txt"
count_file_list=$(find "$rsem_quant_dir" -type f -name "*.expcounts.${strand}.txt" | LC_ALL=C sort | tr '\n' ' ')

CMD="bash ${code_dir}/scripts/join-many.sh $'\t' 0 $count_file_list > ${merged_counts_file}.${sample}.tmp"
(eval "$CMD")
mv -fv "${merged_counts_file}.${sample}.tmp" "$merged_counts_file"


merged_counts_file="${proj_dir}/quant.rsem.TPM.${strand}.txt"
count_file_list=$(find "$rsem_quant_dir" -type f -name "*.TPM.${strand}.txt" | LC_ALL=C sort | tr '\n' ' ')

CMD="bash ${code_dir}/scripts/join-many.sh $'\t' 0 $count_file_list > ${merged_counts_file}.${sample}.tmp"
(eval "$CMD")
mv -fv "${merged_counts_file}.${sample}.tmp" "$merged_counts_file"


merged_counts_file="${proj_dir}/quant.rsem.FPKM.${strand}.txt"
count_file_list=$(find "$rsem_quant_dir" -type f -name "*.FPKM.${strand}.txt" | LC_ALL=C sort | tr '\n' ' ')

CMD="bash ${code_dir}/scripts/join-many.sh $'\t' 0 $count_file_list > ${merged_counts_file}.${sample}.tmp"
(eval "$CMD")
mv -fv "${merged_counts_file}.${sample}.tmp" "$merged_counts_file"


#########################



# end
