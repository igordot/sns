#!/bin/bash


# merge fastqs or create symlinks to originals if no need to merge


# script filename
script_name=$(basename "${BASH_SOURCE[0]}")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 2 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name [project dir] [sample name] \n" >&2
	exit 1
fi

# arguments
proj_dir=$(readlink -f $1)
sample=$2

# check if input exists
if [ ! -d "$proj_dir" ] || [ ! -n "$proj_dir" ] ; then
	echo -e "\n $script_name ERROR: PROJ DIR $proj_dir DOES NOT EXIST \n"
	exit 1
fi

# sample sheet
samples_csv="${proj_dir}/samples.fastq-raw.csv"
if [ ! -s "$samples_csv" ] || [ ! -n "$samples_csv" ] ; then
	echo -e "\n $script_name ERROR: SAMPLES CSV $samples_csv DOES NOT EXIST \n"
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

	# check that R1 file exists and not set to null
	if [ ! -e "$fastq_R1" ] || [ ! -n "$fastq_R1" ] ; then
		echo -e "\n $script_name ERROR: FASTQ R1 $fastq_R1 DOES NOT EXIST \n"
		exit 1
	fi

	# check that R2 file exists if specified
	if [ ! -e "$fastq_R2" ] && [ -n "$fastq_R2" ] ; then
		echo -e "\n $script_name ERROR: FASTQ R2 $fastq_R2 DOES NOT EXIST \n"
		exit 1
	fi

	if [ "$num_files" -gt 1 ] ; then

		# merge multple FASTQs

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


# check if output exists

if [ ! -s "$fastq_R1_clean" ] || [ ! -n "$fastq_R1_clean" ] ; then
	echo -e "\n $script_name ERROR: FASTQ $fastq_R1_clean DOES NOT EXIST \n"
	exit 1
fi


#########################


# count number of reads in merged FASTQs

lines_R1=$(zcat "$fastq_R1_clean" | wc -l)
reads_R1=$(echo "${lines_R1}/4" | bc)

if [ -s "$fastq_R2_clean" ] ; then
	lines_R2=$(zcat "$fastq_R2_clean" | wc -l)
	reads_R2=$(echo "${lines_R2}/4" | bc)
else
	reads_R2="0"
	fastq_R2_clean=""
fi

echo "FASTQ R1: $fastq_R1_clean"
echo "FASTQ R2: $fastq_R2_clean"
echo "READS R1: $reads_R1"
echo "READS R2: $reads_R2"


#########################


# generate summary

echo "#SAMPLE,R1 RAW READS,R2 RAW READS" > "$summary_csv"
echo "${sample},${reads_R1},${reads_R2}" >> "$summary_csv"

sleep 30

# combine all sample summaries
cat ${summary_dir}/*.${segment_name}.csv | LC_ALL=C sort -t ',' -k1,1 | uniq > "${proj_dir}/summary.${segment_name}.csv"


#########################


# exit if FASTQ has very few reads

if [ $reads_R1 -lt 10000 ] ; then
	echo -e "\n $script_name ERROR: FASTQ $fastq_R1_clean IS TOO SHORT \n" >&2
	# delete FASTQs since they are not useable
	rm -fv "$fastq_R1_clean"
	if [ -s "$fastq_R2_clean" ] ; then
		rm -fv "$fastq_R2_clean"
	fi
	exit 1
fi


#########################


# add sample and FASTQ to sample sheet
echo "${sample},${fastq_R1_clean},${fastq_R2_clean}" >> "$samples_csv_clean"

sleep 1

# add again (reduce potential loss if another sample is sorting at the same time)
echo "${sample},${fastq_R1_clean},${fastq_R2_clean}" >> "$samples_csv_clean"

sleep 1

# sort and remove duplicates in place in sample sheet
LC_ALL=C sort -t ',' -k1,1 -u -o "$samples_csv_clean" "$samples_csv_clean"

sleep 30


#########################



# end
