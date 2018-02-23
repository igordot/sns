#!/bin/bash


# merge fastqs or create symlinks to originals if no need to merge


# script filename
script_path="${BASH_SOURCE[0]}"
script_name=$(basename "$script_path")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 2 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir sample_name \n" >&2
	if [ $# -gt 0 ] ; then echo -e "\n ARGS: $* \n" >&2 ; fi
	exit 1
fi

# arguments
proj_dir=$(readlink -f $1)
sample=$2

# check if input exists
if [ ! -d "$proj_dir" ] || [ ! -n "$proj_dir" ] ; then
	echo -e "\n $script_name ERROR: PROJ DIR $proj_dir DOES NOT EXIST \n" >&2
	exit 1
fi

# sample sheet
samples_csv="${proj_dir}/samples.fastq-raw.csv"
if [ ! -s "$samples_csv" ] || [ ! -n "$samples_csv" ] ; then
	echo -e "\n $script_name ERROR: SAMPLES CSV $samples_csv DOES NOT EXIST \n" >&2
	exit 1
fi
samples_csv_clean="${proj_dir}/samples.${segment_name}.csv"

# output summary
summary_dir="${proj_dir}/summary"
mkdir -p "${summary_dir}"
summary_csv="${summary_dir}/${sample}.${segment_name}.csv"

# output files
fastq_clean_dir="${proj_dir}/FASTQ-CLEAN"
mkdir -p "$fastq_clean_dir"
fastq_R1_clean="${fastq_clean_dir}/${sample}_R1.fastq.gz"
fastq_R2_clean="${fastq_clean_dir}/${sample}_R2.fastq.gz"

# check if output already exists
if [ -s "$fastq_R1_clean" ] || [ -s "$fastq_R2_clean" ]; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	echo -e "\n $script_name ADD $sample TO $samples_csv_clean \n" >&2
	if [ -s "$fastq_R2_clean" ] ; then
		echo "${sample},${fastq_R1_clean},${fastq_R2_clean}" >> "$samples_csv_clean"
	else
		echo "${sample},${fastq_R1_clean}," >> "$samples_csv_clean"
	fi
	exit 1
fi

# get number of files (or file pairs) for the sample (to decide if files need to be merged)
num_files=$(grep -c "^${sample}," "${samples_csv}")

# scan for files that belong to the sample
grep "^${sample}," "${samples_csv}" | LC_ALL=C sort | while read -r LINE ; do

	fastq_R1=$(echo "$LINE" | cut -d "," -f 2)
	fastq_R2=$(echo "$LINE" | cut -d "," -f 3)
	fastq_R1=$(readlink -f "$fastq_R1")
	fastq_R2=$(readlink -f "$fastq_R2")

	# check that R1 file is not set to null and exists
	if [ ! -n "$fastq_R1" ] || [ ! -e "$fastq_R1" ] ; then
		echo -e "\n $script_name ERROR: INPUT FASTQ R1 $fastq_R1 DOES NOT EXIST \n" >&2
		exit 1
	fi

	# check that R2 file exists if specified
	if [ -n "$fastq_R2" ] && [ ! -e "$fastq_R2" ] ; then
		echo -e "\n $script_name ERROR: INPUT FASTQ R2 $fastq_R2 DOES NOT EXIST \n" >&2
		exit 1
	fi

	# check R1 file integrity and clear the output file if there is a problem
	if ! gzip --test "$fastq_R1" ; then
		echo -e "\n $script_name ERROR: INPUT FASTQ R1 $fastq_R1 IS CORRUPT \n" >&2
		num_files="0"
		echo "." > "$fastq_R1_clean"
		exit 1
	fi

	# check R2 file integrity if specified and clear the output file if there is a problem
	if [ -n "$fastq_R2" ] ; then
		if ! gzip --test "$fastq_R2" ; then
			echo -e "\n $script_name ERROR: INPUT FASTQ R2 $fastq_R2 IS CORRUPT \n" >&2
			num_files="0"
			echo "." > "$fastq_R2_clean"
			exit 1
		fi
	fi

	if [ "$num_files" -gt 1 ] ; then

		# merge multiple FASTQs

		bash_cmd="cat $fastq_R1 >> $fastq_R1_clean"
		echo "CMD: $bash_cmd"
		eval "$bash_cmd"

		if [ -n "$fastq_R2" ] ; then
			bash_cmd="cat $fastq_R2 >> $fastq_R2_clean"
			echo "CMD: $bash_cmd"
			eval "$bash_cmd"
		fi

	elif [ "$num_files" -eq 1 ] ; then

		# symlink to the original FASTQ if only one

		bash_cmd="ln -s $fastq_R1 $fastq_R1_clean"
		echo "CMD: $bash_cmd"
		eval "$bash_cmd"

		if [ -n "$fastq_R2" ] ; then
			bash_cmd="ln -s $fastq_R2 $fastq_R2_clean"
			echo "CMD: $bash_cmd"
			eval "$bash_cmd"
		fi

	fi

	sleep 30

	date

done


#########################


# count number of reads in clean output FASTQs

reads_R1="0"
reads_R2="0"

if [ -s "$fastq_R1_clean" ] ; then
	lines_R1=$(zcat "$fastq_R1_clean" | wc -l)
	reads_R1=$(echo "${lines_R1}/4" | bc)
fi

if [ -s "$fastq_R2_clean" ] ; then
	lines_R2=$(zcat "$fastq_R2_clean" | wc -l)
	reads_R2=$(echo "${lines_R2}/4" | bc)
else
	fastq_R2_clean=""
fi

echo "FASTQ R1: $fastq_R1_clean"
echo "FASTQ R2: $fastq_R2_clean"
echo "READS R1: $reads_R1"
echo "READS R2: $reads_R2"


#########################


# generate summary
# problematic samples are included so there is a record of them in the summary tables

echo "#SAMPLE,R1 RAW READS,R2 RAW READS" > "$summary_csv"
echo "${sample},${reads_R1},${reads_R2}" >> "$summary_csv"

sleep 30

# combine all sample summaries
cat ${summary_dir}/*.${segment_name}.csv | LC_ALL=C sort -t ',' -k1,1 | uniq > "${proj_dir}/summary.${segment_name}.csv"


#########################


# check FASTQ read numbers

# check if output exists at all
if [ ! -s "$fastq_R1_clean" ] || [ ! -n "$fastq_R1_clean" ] ; then
	echo -e "\n $script_name ERROR: OUTPUT FASTQ $fastq_R1_clean DOES NOT EXIST \n" >&2
	exit 1
fi

# check if there are very few reads
if [ "$reads_R1" -lt 1000 ] ; then
	echo -e "\n $script_name ERROR: OUTPUT FASTQ $fastq_R1_clean IS TOO SHORT \n" >&2
	# delete clean (not original) FASTQs since they are not usable
	rm -fv "$fastq_R1_clean"
	if [ -s "$fastq_R2_clean" ] ; then
		rm -fv "$fastq_R2_clean"
	fi
	exit 1
fi

# check if R1 and R2 have equal number of reads
if [ -s "$fastq_R2_clean" ] ; then
	if [ "$reads_R1" -ne "$reads_R2" ] ; then
		echo -e "\n $script_name ERROR: FASTQ R1 AND R2 HAVE DIFFERENT NUMBER OF READS \n" >&2
		# delete clean (not original) FASTQs since they are not usable
		rm -fv "$fastq_R1_clean"
		rm -fv "$fastq_R2_clean"
		exit 1
	fi
fi


#########################


# add sample and FASTQ to sample sheet
echo "${sample},${fastq_R1_clean},${fastq_R2_clean}" >> "$samples_csv_clean"

sleep 30


#########################



# end
