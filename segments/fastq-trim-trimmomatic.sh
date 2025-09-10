#!/bin/bash


# run Trimmomatic


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
proj_dir=$1
sample=$2
threads=$3
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
fastq_R1_trim_unpaired="${fastq_trim_dir}/${sample}_R1.unpaired.fastq.gz"
fastq_R2_trim_unpaired="${fastq_trim_dir}/${sample}_R2.unpaired.fastq.gz"

trimmomatic_logs_dir="${proj_dir}/logs-trimmomatic"
mkdir -p "$trimmomatic_logs_dir"
trimmomatic_log="${trimmomatic_logs_dir}/${sample}.txt"

# unload all loaded modulefiles
module purge


#########################


# check for output

# delete likely incomplete trimmed FASTQs if unpaired FASTQs exist
if [ -s "$fastq_R1_trim_unpaired" ] ; then
	echo -e "\n $script_name WARNING: POTENTIALLY CORRUPT FASTQ $fastq_R1_trim_unpaired EXISTS \n" >&2
	rm -fv "$fastq_R1_trim"
	rm -fv "$fastq_R2_trim"
	rm -fv "$fastq_R1_trim_unpaired"
	rm -fv "$fastq_R2_trim_unpaired"
fi

# check if output exists already
if [ -s "$fastq_R1_trim" ] && [ -s "$trimmomatic_log" ] ; then

	# check if Trimmomatic reported successfull completion
	if grep -q "Completed successfully" "$trimmomatic_log" ; then
		# exit if output exists already
		echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
		echo -e "\n $script_name ADD $sample TO $samples_csv \n" >&2
		echo "${sample},${fastq_R1_trim},${fastq_R2_trim}" >> "$samples_csv"
		exit 0
	else
		# remove partial output
		echo -e "\n $script_name WARNING: POTENTIALLY INCOMPLETE FASTQ $fastq_R1_trim EXISTS \n" >&2
		rm -fv "$trimmomatic_log"
		rm -fv "$fastq_R1_trim"
		rm -fv "$fastq_R2_trim"
		rm -fv "$fastq_R1_trim_unpaired"
		rm -fv "$fastq_R2_trim_unpaired"
	fi

fi


#########################


# Trimmomatic

module add trimmomatic/0.40

trimmomatic_jar="${TRIMMOMATIC_ROOT}/trimmomatic-0.40.jar"

# check if the trimmomatic jar file is present
if [ ! -s "$trimmomatic_jar" ] ; then
	echo -e "\n $script_name ERROR: FILE $trimmomatic_jar DOES NOT EXIST \n" >&2
	exit 1
fi

if [ -n "$fastq_R2" ] ; then
	run_type_arg="PE"
	files_arg="$fastq_R1 $fastq_R2 $fastq_R1_trim $fastq_R1_trim_unpaired $fastq_R2_trim $fastq_R2_trim_unpaired"
else
	run_type_arg="SE"
	files_arg="$fastq_R1 $fastq_R1_trim"
	fastq_R2_trim=""
	fastq_R2_trim_unpaired=""
fi

echo
echo " * Java bin: $(readlink -f $(which java)) "
echo " * Java version: $(java -version 2>&1 | grep -m 1 'version') "
echo " * Trimmomatic jar: $trimmomatic_jar "
echo " * Trimmomatic version: $(java -jar $trimmomatic_jar -version) "
echo " * run type: $run_type_arg "
echo " * FASTQ R1: $fastq_R1 "
echo " * FASTQ R2: $fastq_R2 "
echo " * FASTQ R1 trimmed: $fastq_R1_trim "
echo " * FASTQ R2 trimmed: $fastq_R2_trim "
echo " * FASTQ R1 trimmed unpaired: $fastq_R1_trim_unpaired "
echo " * FASTQ R2 trimmed unpaired: $fastq_R2_trim_unpaired "
echo

bash_cmd="
java -Xms32G -Xmx32G -jar $trimmomatic_jar \
$run_type_arg \
-threads $threads \
$files_arg \
ILLUMINACLIP:/gpfs/data/igorlab/ref/contaminants/trimmomatic.040.fa:2:30:10:2:True \
TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:35 \
2> $trimmomatic_log
"
echo "CMD: $bash_cmd"
eval "$bash_cmd"

sleep 5

# delete unpaired files
rm -fv "$fastq_R1_trim_unpaired"
rm -fv "$fastq_R2_trim_unpaired"


#########################


# check that output generated

if [ ! -s "$fastq_R1_trim" ] ; then
	echo -e "\n $script_name ERROR: FASTQ $fastq_R1_trim NOT GENERATED \n" >&2
	exit 1
fi

if grep -q -i "error" "$trimmomatic_log" ; then
	echo -e "\n $script_name ERROR: LOG $trimmomatic_log CONTAINS ERROR \n" >&2
	grep -i "error" "$trimmomatic_log"
	rm -fv "$fastq_R1_trim"
	rm -fv "$fastq_R2_trim"
	exit 1
fi


#########################


# generate summary

# extract the relevant part of the log

if [ -n "$fastq_R2" ] ; then
	total_reads_col="4"
	surviving_reads_col="7"
else
	total_reads_col="3"
	surviving_reads_col="5"
fi

reads_input=$(cat $trimmomatic_log | grep "Input Read" | cut -d ' ' -f $total_reads_col)
reads_surviving=$(cat $trimmomatic_log | grep "Input Read" | cut -d ' ' -f $surviving_reads_col)

reads_surviving_pct=$(echo "(${reads_surviving}/${reads_input})*100" | bc -l | cut -c 1-4)
reads_surviving_pct="${reads_surviving_pct}%"

# header for summary file
echo "#SAMPLE,TRIM INPUT READS,TRIM SURVIVING READS,TRIM SURVIVING READS %" > "$summary_csv"

# summarize log file
echo "${sample},${reads_input},${reads_surviving},${reads_surviving_pct}" >> "$summary_csv"

sleep 5

# combine all sample summaries
cat ${summary_dir}/*.${segment_name}.csv | LC_ALL=C sort -t ',' -k1,1 | uniq > "${proj_dir}/summary.${segment_name}.csv"


#########################


# add sample and FASTQ to sample sheet
echo "${sample},${fastq_R1_trim},${fastq_R2_trim}" >> "$samples_csv"

sleep 5


#########################



# end
