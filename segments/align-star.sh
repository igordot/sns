#!/bin/bash


# run STAR


# script filename
script_path="${BASH_SOURCE[0]}"
script_name=$(basename "$script_path")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ $# -lt 4 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir sample_name num_threads FASTQ_R1 [FASTQ_R2] \n" >&2
	if [ $# -gt 0 ] ; then echo -e "\n ARGS: $* \n" >&2 ; fi
	exit 1
fi

# arguments
proj_dir=$(readlink -f "$1")
sample=$2
threads=$3
fastq_R1=$4
fastq_R2=$5


#########################


# settings and files

summary_dir="${proj_dir}/summary"
mkdir -p "$summary_dir"
summary_csv="${summary_dir}/${sample}.${segment_name}.csv"

samples_csv="${proj_dir}/samples.${segment_name}.csv"

star_bam_dir="${proj_dir}/BAM-STAR"
mkdir -p "$star_bam_dir"
bam="${star_bam_dir}/${sample}.bam"
bai="${bam}.bai"

star_quant_dir="${proj_dir}/quant-STAR"
mkdir -p "$star_quant_dir"
counts_txt="${star_quant_dir}/${sample}.txt"

star_logs_dir="${proj_dir}/logs-${segment_name}"
mkdir -p "$star_logs_dir"
star_prefix="${star_logs_dir}/${sample}_"

# unload all loaded modulefiles
module purge


#########################


# check for output

# skip if final BAM and BAI exist
if [ -s "$bam" ] && [ -s "$bai" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	echo -e "\n $script_name ADD $sample TO $samples_csv \n" >&2
	echo "${sample},${bam}" >> "$samples_csv"
	exit 0
fi

# delete BAM (likely incomplete since the corresponding BAI was not generated)
if [ -s "$bam" ] ; then
	echo -e "\n $script_name WARNING: potentially corrupt BAM $bam exists \n" >&2
	rm -fv "$bam"
fi

# skip if run is in progress
if [ -s "${star_prefix}.bam" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 0
fi


#########################


# check that inputs exist

if [ ! -d "$proj_dir" ] ; then
	echo -e "\n $script_name ERROR: proj dir $proj_dir does not exist \n" >&2
	exit 1
fi

if [ ! -s "$fastq_R1" ] ; then
	echo -e "\n $script_name ERROR: FASTQ $fastq_R1 does not exist \n" >&2
	exit 1
fi

code_dir=$(dirname $(dirname "$script_path"))

ref_star=$(bash "${code_dir}/scripts/get-set-setting.sh" "${proj_dir}/settings.txt" REF-STAR)

if [ ! -s "${ref_star}/SA" ] ; then
	echo -e "\n $script_name ERROR: STAR reference $ref_star does not exist \n" >&2
	exit 1
fi


#########################


# STAR

# star/2.7.3a loads samtools/1.9
# reference genome version 2.7.1a (updated in 2.7.4a)
module add star/2.7.3a

echo
echo " * STAR: $(readlink -f $(which STAR)) "
echo " * samtools: $(readlink -f $(which samtools)) "
echo " * samtools version: $(samtools --version | head -1) "
echo " * STAR ref: $ref_star "
echo " * FASTQ R1: $fastq_R1 "
echo " * FASTQ R2: $fastq_R2 "
echo " * BAM: $bam "
echo

# change dir because STAR and samtools may generate temp or log files in working dir
cd "$star_logs_dir"

# --outFilterType BySJout - ENCODE standard option
# --outSAMstrandField intronMotif - add XS strand tags to spliced reads
# --outSAMattrRGline - first word contains the read group identifier and must start with "ID:"
# --outSAMmapqUnique - int: 0 to 255: the MAPQ value for unique mappers

# relevant errors:
# "HTseq cannot deal with jM:B:c,-1 and jI:B:i,-1 SAM attributes which are output if you use (non-default) --outSAMmode Full"
# according to Picard, STAR-sorted BAM "is not coordinate sorted", so using samtools for sorting

bash_cmd="
STAR \
--runThreadN ${threads} \
--genomeDir ${ref_star} \
--genomeLoad NoSharedMemory \
--outFilterMismatchNoverLmax 0.2 \
--outFilterMultimapNmax 1 \
--outFilterType BySJout \
--outSAMstrandField intronMotif \
--outSAMattributes NH HI AS nM NM MD XS \
--outSAMattrRGline ID:${sample} SM:${sample} LB:${sample} PL:ILLUMINA \
--outSAMmapqUnique 60 \
--twopassMode Basic \
--readFilesCommand zcat \
--readFilesIn ${fastq_R1} ${fastq_R2} \
--outFileNamePrefix ${star_prefix} \
--quantMode GeneCounts \
--outSAMtype BAM Unsorted \
--outStd BAM_Unsorted | \
samtools sort -m 32G -T ${sample}.samtools -o ${bam} -
"
echo "CMD: $bash_cmd"
eval "$bash_cmd"

sleep 5


#########################


# delete tmp directories (_STARtmp should be empty at this point)
rm -rfv ${star_prefix}_STAR*


#########################


# check that output generated

# check if BAM file is present
if [ ! -s "$bam" ] ; then
	echo -e "\n $script_name ERROR: BAM $bam not generated \n" >&2
	exit 1
fi

# check if gene counts file is present
if [ ! -s "${star_prefix}ReadsPerGene.out.tab" ] ; then
	echo -e "\n $script_name ERROR: counts ${star_prefix}ReadsPerGene.out.tab not generated \n" >&2
	# delete BAM and BAI since something went wrong and they might be corrupted
	rm -fv "$bam"
	exit 1
fi


#########################


# index BAM

bash_cmd="samtools index $bam"
echo "CMD: $bash_cmd"
eval "$bash_cmd"

sleep 5


#########################


# check that output generated

# check if BAM index is present (generated only if BAM is valid)
if [ ! -s "$bai" ] ; then
	echo -e "\n $script_name ERROR: BAM index $bai not generated \n" >&2
	# delete BAM since something went wrong and it might be corrupted
	rm -fv "$bam"
	exit 1
fi


#########################


# STAR counts
# not using the counts downstream since the file is hard to filter and uses gene ids instead of gene names

# STAR outputs read counts per gene into ReadsPerGene.out.tab file with 4 columns which correspond to strandedness:
# column 1: gene ID
# column 2: counts for unstranded RNA-seq
# column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
# column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)

# move counts file to a separate bam-only directory
CMD="mv -v ${star_prefix}ReadsPerGene.out.tab $counts_txt"
echo "CMD: $CMD"
eval "$CMD"

# strand:
# fwd | transcript             | cufflinks "fr-secondstrand" | htseq "yes"     | picard "FIRST_READ"
# rev | rev comp of transcript | cufflinks "fr-firststrand"  | htseq "reverse" | picard "SECOND_READ"

# top 1000 genes (useful to determine strandness)
echo -e "#gene_id,unstr,fwd,rev" > "${star_quant_dir}/${sample}.top1000.csv"
cat "$counts_txt" | LC_ALL=C sort -k2,2nr | head -n 1000 | tr '\t' ',' >> "${star_quant_dir}/${sample}.top1000.csv"


#########################


# determine strand

# get total counts for each strand
counts_unstr=$(cat $counts_txt | grep -v 'N_' | awk -F $'\t' '{sum+=$2} END {print sum}')
echo "counts unstr: $counts_unstr"
counts_fwd=$(cat $counts_txt | grep -v 'N_' | awk -F $'\t' '{sum+=$3} END {print sum}')
echo "counts fwd: $counts_fwd"
counts_rev=$(cat $counts_txt | grep -v 'N_' | awk -F $'\t' '{sum+=$4} END {print sum}')
echo "counts rev: $counts_rev"

lib_strand="unstr"

if [ "$(echo "${counts_fwd}/${counts_rev}" | bc)" -gt 5 ] ; then
	lib_strand="fwd"
fi

if [ "$(echo "${counts_rev}/${counts_fwd}" | bc)" -gt 5 ] ; then
	lib_strand="rev"
fi

# set experiment strand if the file meets some quality standards (strand not set if previously set)
if [ "$counts_unstr" -gt 10000 ] ; then
	exp_strand=$(bash "${code_dir}/scripts/get-set-setting.sh" "${proj_dir}/settings.txt" EXP-STRAND "$lib_strand");
fi

# generate an error for low quality files
if [ "$counts_unstr" -lt 1000 ] || [ "$counts_fwd" -lt 10 ] || [ "$counts_rev" -lt 10 ] ; then
	echo -e "\n $script_name ERROR: low counts \n" >&2
fi

echo "sample strand: $lib_strand"
echo "experiment strand: $exp_strand"


#########################


# generate alignment summary

# header for summary file
echo "#SAMPLE,INPUT READS,UNIQUELY MAPPED,MULTI-MAPPED,UNIQUELY MAPPED %,MULTI-MAPPED %" > "$summary_csv"

# print the relevant numbers from log file
star_log_final="${star_prefix}Log.final.out"
paste -d ',' \
<(echo "$sample") \
<(cat "$star_log_final" | grep "Number of input reads"                   | head -1 | tr -d "[:blank:]" | cut -d "|" -f 2) \
<(cat "$star_log_final" | grep "Uniquely mapped reads number"            | head -1 | tr -d "[:blank:]" | cut -d "|" -f 2) \
<(cat "$star_log_final" | grep "Number of reads mapped to too many loci" | head -1 | tr -d "[:blank:]" | cut -d "|" -f 2) \
<(cat "$star_log_final" | grep "Uniquely mapped reads %"                 | head -1 | tr -d "[:blank:]" | cut -d "|" -f 2) \
<(cat "$star_log_final" | grep "% of reads mapped to too many loci"      | head -1 | tr -d "[:blank:]" | cut -d "|" -f 2) \
>> "$summary_csv"

sleep 5

# combine all sample summaries
cat ${summary_dir}/*.${segment_name}.csv | LC_ALL=C sort -t ',' -k1,1 | uniq > "${proj_dir}/summary.${segment_name}.csv"


#########################


# add sample and BAM to sample sheet
echo "${sample},${bam}" >> "$samples_csv"

sleep 5


#########################



# end
