#!/bin/bash


# run featureCounts


# script filename
script_name=$(basename "${BASH_SOURCE[0]}")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 6 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir sample_name threads BAM PE/SE strand \n" >&2
	exit 1
fi

# arguments
proj_dir=$1
sample=$2
threads=$3
bam=$4
run_type=$5
strand=$6


#########################


# check that inputs exist

if [ ! -d "$proj_dir" ] || [ ! "$proj_dir" ] ; then
	echo -e "\n $script_name ERROR: DIR $proj_dir DOES NOT EXIST \n" >&2
	exit 1
fi

if [ ! -s "$bam" ] || [ ! "$bam" ] ; then
	echo -e "\n $script_name ERROR: bam $bam DOES NOT EXIST \n" >&2
	exit 1
fi

code_dir=$(dirname "$(dirname "${BASH_SOURCE[0]}")")
gtf=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" REF-GTF);

if [ ! -s "$gtf" ] || [ ! "$gtf" ] ; then
	echo -e "\n $script_name ERROR: GTF $gtf DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# settings and files

segment_name="${segment_name}-${strand}"

summary_dir="${proj_dir}/summary"
mkdir -p "$summary_dir"
summary_csv="${summary_dir}/${sample}.${segment_name}.csv"

fc_quant_dir=${proj_dir}/quant-featureCounts
mkdir -p "$fc_quant_dir"
counts_clean=${fc_quant_dir}/${sample}.counts.${strand}.txt

fc_logs_dir="${proj_dir}/logs-featureCounts"
mkdir -p "$fc_logs_dir"
counts_raw="${fc_logs_dir}/${sample}.${strand}.txt"
counts_summary_raw="${fc_logs_dir}/${sample}.${strand}.txt.summary"

# run type flag
if [ "$run_type" == "SE" ] ; then
	run_type_flag=""
elif [ "$run_type" == "PE" ] ; then
	run_type_flag="-p"
else
	echo -e "\n $script_name ERROR: incorrect run type $run_type selected \n" >&2
	exit 1
fi

# strand flag
if [ "$strand" == "unstr" ] ; then
	strand_flag="0"
elif [ "$strand" == "fwd" ] ; then
	strand_flag="1"
elif [ "$strand" == "rev" ] ; then
	strand_flag="2"
else
	echo -e "\n $script_name ERROR: incorrect strand $strand selected \n" >&2
	exit 1
fi


#########################


# exit if output exits already

if [ -s "$counts_clean" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 1
fi


#########################


# featureCounts

module load subread/1.4.6-p3

# featureCounts generates temp files in the current directory
cd "$fc_logs_dir"

echo " * featureCounts: $(readlink -f $(which featureCounts)) "
echo " * run type: $run_type "
echo " * BAM: $bam "
echo " * GTF: $gtf "
echo " * counts: $counts_clean "

# strand:
# fwd | transcript             | cufflinks "fr-secondstrand" | htseq "yes"     | picard "FIRST_READ"
# rev | rev comp of transcript | cufflinks "fr-firststrand"  | htseq "reverse" | picard "SECOND_READ"

# Summarize a single-end read dataset using 5 threads:
# featureCounts -T 5 -t exon -g gene_id -a annotation.gtf -o counts.txt mapping_results_SE.sam
# Summarize a bam format dataset:
# featureCounts -t exon -g gene_id -a annotation.gtf -o counts.txt mapping_results_SE.bam
# Summarize multiple datasets at the same time:
# featureCounts -t exon -g gene_id -a annotation.gtf -o counts.txt library1.bam library2.bam library3.bam
# Perform strand-specific read counting (use '-s 2' if reversely stranded):
# featureCounts -s 1 -t exon -g gene_id -a annotation.gtf -o counts.txt mapping_results_SE.bam
# Summarize paired-end reads and count fragments (instead of reads):
# featureCounts -p -t exon -g gene_id -a annotation.gtf -o counts.txt mapping_results_PE.bam
# Summarize multiple paired-end datasets:
# featureCounts -p -t exon -g gene_id -a annotation.gtf -o counts.txt library1.bam library2.bam library3.bam
# Count the fragments that have fragment length between 50bp and 600bp only:
# featureCounts -p -P -d 50 -D 600 -t exon -g gene_id -a annotation.gtf -o counts.txt mapping_results_PE.bam
# Count those fragments that have both ends mapped only:
# featureCounts -p -B -t exon -g gene_id -a annotation.gtf -o counts.txt mapping_results_PE.bam
# Exclude chimeric fragments from fragment counting:
# featureCounts -p -C -t exon -g gene_id -a annotation.gtf -o counts.txt mapping_results_PE.bam

CMD="
featureCounts \
-T $threads \
$run_type_flag \
-g gene_name \
-s $strand_flag \
-a $gtf \
-o $counts_raw \
$bam
"
echo "CMD: $CMD"
eval "$CMD"


#########################


# check that output generated

if [ ! -s "$counts_raw" ] ; then
	echo -e "\n $script_name ERROR: COUNTS $counts_raw NOT GENERATED \n" >&2
	exit 1
fi


#########################


# clean up the full output

# 1: Geneid
# 2: Chr
# 3: Start
# 4: End
# 5: strand
# 6: Length
# 7: count

# gene counts table
echo -e "#GENE\t${sample}" > "$counts_clean"
cat $counts_raw | grep -v '#' | grep -v 'Geneid' | cut -f 1,7 | LC_ALL=C sort -k1,1 >> "$counts_clean"

# gene info table
gene_info_file="${proj_dir}/genes.featurecounts.txt"
echo -e "#GENE\tCHR\tSTART\tEND\tstrand\tLENGTH" > "$gene_info_file"
cat $counts_raw | grep -v '#' | grep -v 'Geneid' | cut -f 1-6 | LC_ALL=C sort -k1,1 >> "$gene_info_file"

# delete the full original counts table
rm -fv $counts_raw


#########################


# generate summary

counts_assigned=$(cat $counts_summary_raw | grep -m 1 "^Assigned" | cut -f 2)
echo "COUNTS ASSIGNED: $counts_assigned"

# header for summary file
echo "#SAMPLE,ASSIGNED COUNTS ${strand}" > "$summary_csv"

# summarize log file
echo "${sample},${counts_assigned}" >> "$summary_csv"

sleep 30

# combine all sample summaries
cat ${summary_dir}/*.${segment_name}.csv | LC_ALL=C sort -t ',' -k1,1 | uniq > "${proj_dir}/summary.${segment_name}.csv"


#########################


# generate counts matrix for all samples

# merged counts filename
merged_counts_file="${proj_dir}/quant.featurecounts.counts.${strand}.txt"

# list of all counts files
count_file_list=$(find "$fc_quant_dir" -type f -name "*.counts.${strand}.txt" | LC_ALL=C sort | tr '\n' ' ')

# merged counts filename with sample included to avoid conflicts during join process
CMD="bash ${code_dir}/scripts/join-many.sh $'\t' 0 $count_file_list > ${merged_counts_file}.${sample}.tmp"
(eval "$CMD")

# rewrite any existing merged counts table with current table
mv -fv "${merged_counts_file}.${sample}.tmp" "$merged_counts_file"


#########################



# end
