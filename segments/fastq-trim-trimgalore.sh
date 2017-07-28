#!/bin/bash


# run Trim Galore (for MspI digested RRBS samples)


# script filename
script_name=$(basename "${BASH_SOURCE[0]}")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ $# -lt 4 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir sample_name run_type FASTQ_R1 [FASTQ_R2] \n" >&2
	exit 1
fi

# arguments
proj_dir=$1
sample=$2
run_type=$3
fastq_R1=$4
fastq_R2=$5


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


#########################


# settings and files

summary_dir="${proj_dir}/summary"
mkdir -p "$summary_dir"
summary_csv="${summary_dir}/${sample}.${segment_name}.csv"

samples_csv="${proj_dir}/samples.${segment_name}.csv"

fastq_trim_dir="${proj_dir}/FASTQ-TRIMMED"
mkdir -p "$fastq_trim_dir"
fastq_R1_trim="${fastq_trim_dir}/${sample}_R1.trim.fastq.gz"
fastq_R2_trim="${fastq_trim_dir}/${sample}_R2.trim.fastq.gz"

trim_galore_logs_dir="${proj_dir}/logs-trimgalore"
mkdir -p "$trim_galore_logs_dir"

# adjust for single/paired end
if [ -n "$fastq_R2" ] ; then
	fastq_R1_trim_original=${trim_galore_logs_dir}/$(basename $fastq_R1 | sed 's/.fastq.gz/_val_1.fq.gz/')
	fastq_R2_trim_original=${trim_galore_logs_dir}/$(basename $fastq_R2 | sed 's/.fastq.gz/_val_2.fq.gz/')
	trim_galore_log=${trim_galore_logs_dir}/$(basename $fastq_R2)_trimming_report.txt
else
	fastq_R1_trim_original=${trim_galore_logs_dir}/$(basename $fastq_R1 | sed 's/.fastq.gz/_trimmed.fq.gz/')
	trim_galore_log=${trim_galore_logs_dir}/$(basename $fastq_R1)_trimming_report.txt
	fastq_R2_trim=""
fi


#########################


# exit if output exists already

if [ -s "$fastq_R1_trim" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	echo "${sample},${fastq_R1_trim},${fastq_R2_trim}" >> "$samples_csv"
	exit 1
fi


#########################


# trim galore

module load trim-galore/0.4.1

echo " * trim_galore: $(readlink -f $(which trim_galore)) "
echo " * trim_galore version: $(trim_galore --version | grep 'version' | tr -s '[:blank:]') "
echo " * cutadapt: $(readlink -f $(which cutadapt)) "
echo " * cutadapt version: $(cutadapt --version) "
echo " * FASTQ R1: $fastq_R1 "
echo " * FASTQ R2: $fastq_R2 "
echo " * FASTQ R1 TRIMMED: $fastq_R1_trim "
echo " * FASTQ R2 TRIMMED: $fastq_R2_trim "

# --rrbs Specifies that the input file was an MspI digested RRBS sample (recognition site: CCGG)

# optional fastqc
# --fastqc --fastqc_args ' --threads $THREADS --outdir $FASTQC_DIR ' \

# check for run type
if [ "$run_type" == "rrbs" ] ; then
	trim_galore_flags="--length 30 --rrbs"
else
	echo -e "\n $script_name unknown run type $run_type \n" >&2
	exit 1
fi

# adjust for single/paired end
if [ -n "$fastq_R2" ] ; then
	trim_galore_flags="$trim_galore_flags --paired $fastq_R1 $fastq_R2"
else
	trim_galore_flags="$trim_galore_flags $fastq_R1"
fi

# trim_galore
bash_cmd="
trim_galore \
--output_dir $trim_galore_logs_dir \
$trim_galore_flags
"
echo "CMD: $bash_cmd"
eval "$bash_cmd"

sleep 30

# clean up file names

mv -fv "$fastq_R1_trim_original" "$fastq_R1_trim"

if [ -n "$fastq_R2" ] ; then
	mv -fv "$fastq_R2_trim_original" "$fastq_R2_trim"
fi

sleep 30


#########################


# check that output generated

if [ ! -s "$fastq_R1_trim" ] ; then
	echo -e "\n $script_name ERROR: FASTQ $fastq_R1_trim NOT GENERATED \n" >&2
	exit 1
fi


#########################


# generate summary

# extract the relevant part of the log
reads_input=$(cat "$trim_galore_log" | grep "sequences processed in total" | tail -1 | cut -d " " -f 1)

# single-end and paired-end trim logs have different formatting, so just manually counting number of lines
reads_surviving=$(zcat $fastq_R1_trim | wc -l)
reads_surviving=$(echo "${reads_surviving}/4" | bc)

reads_surviving_pct=$(echo "(${reads_surviving}/${reads_input})*100" | bc -l | cut -c 1-4)
reads_surviving_pct="${reads_surviving_pct}%"

# header for summary file
echo "#SAMPLE,TRIM INPUT READS,TRIM SURVIVING READS,TRIM SURVIVING READS %" > "$summary_csv"

# summarize log file
echo "${sample},${reads_input},${reads_surviving},${reads_surviving_pct}" >> "$summary_csv"

sleep 30

# combine all sample summaries
cat ${summary_dir}/*.${segment_name}.csv | LC_ALL=C sort -t ',' -k1,1 | uniq > "${proj_dir}/summary.${segment_name}.csv"


#########################


# add sample and FASTQ to sample sheet
echo "${sample},${fastq_R1_trim},${fastq_R2_trim}" >> "$samples_csv"

sleep 30


#########################



# end
