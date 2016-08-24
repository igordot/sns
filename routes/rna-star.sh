#!/bin/bash


##
## RNA-seq using STAR aligner
##


# script filename
script_name=$(basename "${BASH_SOURCE[0]}")
route_name=${script_name/%.sh/}
echo -e "\n ========== ROUTE: $route_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 2 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name [project dir] [sample name] \n" >&2
	exit 1
fi

# standard route arguments
proj_dir=$(readlink -f "$1")
sample=$2

# additional settings
threads=$NSLOTS
code_dir=$(dirname "$(dirname "${BASH_SOURCE[0]}")")
qsub_dir="${proj_dir}/logs-qsub"

# display settings
echo " * proj_dir: $proj_dir "
echo " * sample: $sample "
echo " * code_dir: $code_dir "
echo " * qsub_dir: $qsub_dir "
echo " * threads: $threads "


#########################


# delete empty qsub .po files
rm -f ${qsub_dir}/sns.*.po*


#########################


# segments

# rename and/or merge raw input FASTQs
segment_fastq_clean="fastq-fastq-clean"
fastq_R1=$(grep -s -m 1 "^${sample}," "${proj_dir}/samples.${segment_fastq_clean}.csv" | cut -d ',' -f 2)
fastq_R2=$(grep -s -m 1 "^${sample}," "${proj_dir}/samples.${segment_fastq_clean}.csv" | cut -d ',' -f 3)
if [ -z "$fastq_R1" ] ; then
	bash_cmd="bash ${code_dir}/segments/${segment_fastq_clean}.sh $proj_dir $sample"
	($bash_cmd)
	fastq_R1=$(grep -m 1 "^${sample}," "${proj_dir}/samples.${segment_fastq_clean}.csv" | cut -d ',' -f 2)
	fastq_R2=$(grep -m 1 "^${sample}," "${proj_dir}/samples.${segment_fastq_clean}.csv" | cut -d ',' -f 3)
fi

# fastq_screen
bash_cmd="bash ${code_dir}/segments/fastq-qc-fastqscreen.sh $proj_dir $sample $fastq_R1"
($bash_cmd)

# run STAR
segment_align="fastq-bam-star"
bam_star=$(grep -s -m 1 "^${sample}," "${proj_dir}/samples.bam-star.csv" | cut -d ',' -f 2)
if [ -z "$bam_star" ] ; then
	bash_cmd="bash ${code_dir}/segments/${segment_align}.sh $proj_dir $sample $threads $fastq_R1 $fastq_R2"
	($bash_cmd)
	bam_star=$(grep -m 1 "^${sample}," "${proj_dir}/samples.${segment_align}.csv" | cut -d ',' -f 2)
fi

# generate BigWig (deeptools)
segment_bigwig_deeptools="bam-bigwig-deeptools"
bash_cmd="bash ${code_dir}/segments/${segment_bigwig_deeptools}.sh $proj_dir $sample 4 $bam_star"
qsub_cmd="qsub -N sns.${segment_bigwig_deeptools}.${sample} -M ${USER}@nyumc.org -m a -j y -cwd -pe threaded 4 -b y ${bash_cmd}"
$qsub_cmd

# generate BigWig (bedtools)
segment_bigwig_bedtools="bam-bigwig-bedtools"
bash_cmd="bash ${code_dir}/segments/${segment_bigwig_bedtools}.sh $proj_dir $sample $bam_star"
qsub_cmd="qsub -N sns.${segment_bigwig_bedtools}.${sample} -M ${USER}@nyumc.org -m a -j y -cwd -b y ${bash_cmd}"
$qsub_cmd

# Picard CollectRnaSeqMetrics
segment_qc_picard="bam-qc-picard-rnaseqmetrics"
bash_cmd="bash ${code_dir}/segments/${segment_qc_picard}.sh $proj_dir $sample $bam_star"
($bash_cmd)

# determine run type for featurecounts
if [ -n "$fastq_R2" ] ; then
	run_type="PE"
else
	run_type="SE"
fi

# determine strand
exp_strand=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" EXP-STRAND);

# generate counts
segment_quant="bam-quant-featurecounts"
bash_cmd="bash ${code_dir}/segments/${segment_quant}.sh $proj_dir $sample $threads $bam_star $run_type $exp_strand"
($bash_cmd)

# generate unstranded counts just in case
if [ "$exp_strand" != "unstr" ] ; then
	bash_cmd="bash ${code_dir}/segments/${segment_quant}.sh $proj_dir $sample $threads $bam_star $run_type unstr"
	($bash_cmd)
fi


#########################


# combine summary from each step

sleep 30

summary_csv="${proj_dir}/summary-combined.rna-star.csv"

bash_cmd="
bash ${code_dir}/scripts/join-many.sh , X \
${proj_dir}/summary.${segment_fastq_clean}.csv \
${proj_dir}/summary.${segment_align}.csv \
${proj_dir}/summary.${segment_quant}-unstr.csv \
${proj_dir}/summary.${segment_quant}-${exp_strand}.csv \
${proj_dir}/summary.${segment_qc_picard}.csv \
> $summary_csv
"
(eval $bash_cmd)


#########################


# generate groups sample sheet template

samples_groups_csv="${proj_dir}/samples.groups.csv"

if [ ! -s "$samples_groups_csv" ] ; then
	echo "#SAMPLE,group" > $samples_groups_csv
	sed 's/\,.*/,NA/g' samples.fastq-raw.csv | LC_ALL=C sort -u >> $samples_groups_csv
fi


#########################


# delete empty qsub .po files
rm -f ${qsub_dir}/sns.*.po*


#########################


date



# end
