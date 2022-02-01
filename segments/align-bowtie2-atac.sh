#!/bin/bash


# run Bowtie2 with ATAC-seq specific parameters


# script filename
script_path="${BASH_SOURCE[0]}"
script_name=$(basename "$script_path")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 5 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir sample_name num_threads FASTQ_R1 FASTQ_R2 \n" >&2
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


# check that inputs exist

if [ ! -d "$proj_dir" ] ; then
	echo -e "\n $script_name ERROR: DIR $proj_dir DOES NOT EXIST \n" >&2
	exit 1
fi

if [ ! -s "$fastq_R1" ] ; then
	echo -e "\n $script_name ERROR: FASTQ $fastq_R1 DOES NOT EXIST \n" >&2
	exit 1
fi

if [ ! -s "$fastq_R2" ] ; then
	echo -e "\n $script_name ERROR: FASTQ $fastq_R2 DOES NOT EXIST \n" >&2
	exit 1
fi

code_dir=$(dirname $(dirname "$script_path"))

ref_bowtie2=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" REF-BOWTIE2);

if [ ! -s "${ref_bowtie2}.fa" ] ; then
	echo -e "\n $script_name ERROR: REF $ref_bowtie2 DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# settings and files

summary_dir="${proj_dir}/summary"
mkdir -p "$summary_dir"
summary_csv="${summary_dir}/${sample}.${segment_name}.csv"

samples_csv="${proj_dir}/samples.${segment_name}.csv"

bam_dir="${proj_dir}/BAM"
mkdir -p "$bam_dir"
bam="${bam_dir}/${sample}.bam"
bai="${bam}.bai"

logs_dir="${proj_dir}/logs-${segment_name}"
mkdir -p "$logs_dir"
unfiltered_sam="${logs_dir}/${sample}.unfiltered.sam"
bowtie2_txt="${logs_dir}/${sample}.bowtie2.txt"
flagstat_txt="${logs_dir}/${sample}.flagstat.txt"

# unload all loaded modulefiles
module purge
module add default-environment


#########################


# check for output

# skip if final BAM and BAI exist
if [ -s "$bam" ] && [ -s "$bai" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	echo -e "\n $script_name ADD $sample TO $samples_csv \n" >&2
	echo "${sample},${bam}" >> "$samples_csv"
	exit 1
fi

# delete BAM (likely incomplete since the corresponding BAI was not generated)
if [ -s "$bam" ] ; then
	echo -e "\n $script_name WARNING: POTENTIALLY CORRUPT BAM $bam EXISTS \n" >&2
	rm -fv "$bam"
fi

# delete SAM (likely incomplete since the BAM was not generated)
if [ -s "$unfiltered_sam" ] ; then
	echo -e "\n $script_name WARNING: SAM $unfiltered_sam EXISTS \n" >&2
	rm -fv "$unfiltered_sam"
fi


#########################


# Bowtie2

# step 1: align with Bowtie2
# step 2: convert SAM to BAM and remove low quality reads
# step 3: sort BAM
# to do: add Picard AddOrReplaceReadGroups for extra compatibility

module add bowtie2/2.3.4.1
module add samtools/1.9
module add sambamba/0.6.8

sambamba_bin="sambamba-0.6.8"

echo
echo " * bowtie2: $(readlink -f $(which bowtie2)) "
echo " * bowtie2 version: $(bowtie2 --version 2>&1 | head -1) "
echo " * sambamba: $(readlink -f $(which $sambamba_bin)) "
echo " * sambamba version: $($sambamba_bin 2>&1 | grep -m 1 'sambamba') "
echo " * samtools: $(readlink -f $(which samtools)) "
echo " * samtools version: $(samtools --version | head -1) "
echo " * bowtie2 ref: $ref_bowtie2 "
echo " * FASTQ R1: $fastq_R1 "
echo " * FASTQ R2: $fastq_R2 "
echo " * BAM: $bam "
echo

# --no-mixed
#  by default, bowtie2 tries to align each mate if it cannot find a concordant or discordant alignment for a pair
# --no-discordant
#  a discordant alignment does not satisfy the paired-end constraints (--fr/--rf/--ff, -I, -X)
# --dovetail
#  if one mate alignment extends past the beginning of the other, consider that to be concordant
# --soft-clipped-unmapped-tlen
#  consider soft-clipped bases unmapped when calculating TLEN

bash_cmd="
bowtie2 \
--local \
--minins 25 \
--maxins 2000 \
--no-mixed \
--no-discordant \
--dovetail \
--soft-clipped-unmapped-tlen \
--threads $threads \
-x $ref_bowtie2 \
-1 $fastq_R1 -2 $fastq_R2 \
-S $unfiltered_sam \
2> $bowtie2_txt \
"
echo "CMD: $bash_cmd"
eval "$bash_cmd"


#########################


# check that output generated

# check that file size above 10 kb
sam_size=$(du -s --apparent-size $unfiltered_sam | cut -f 1)
if [ $sam_size -lt 10 ] ; then
	echo -e "\n $script_name ERROR: unfiltered SAM $unfiltered_sam too small \n" >&2
	exit 1
fi


#########################


# get pre-filtered alignment stats

reads_input=$(samtools view -c "$unfiltered_sam")
echo "READS INPUT: $reads_input"

reads_mapped=$(samtools view -c -F 4 "$unfiltered_sam")
echo "READS MAPPED: $reads_mapped"

reads_chrM=$(samtools view "$unfiltered_sam" | cut -f 3 | grep -c "chrM")
echo "READS CHR M: $reads_chrM"


#########################


# filter and sort BAM

bash_cmd="
cat $unfiltered_sam \
| \
$sambamba_bin view \
--sam-input \
--nthreads $threads \
--filter \"mapping_quality >= 30 and ref_name != 'chrM'\" \
--format bam \
--compression-level 0 \
/dev/stdin \
| \
$sambamba_bin sort \
--nthreads $threads \
--memory-limit 8GB \
--out $bam \
/dev/stdin
"
echo "CMD: $bash_cmd"
eval "$bash_cmd"

rm -fv "$unfiltered_sam"


#########################


# check that output generated

# check if BAM file is present
if [ ! -s "$bam" ] ; then
	echo -e "\n $script_name ERROR: BAM $bam NOT GENERATED \n" >&2
	exit 1
fi

# check if BAM index is present (generated only if BAM is valid)
if [ ! -s "$bai" ] ; then
	echo -e "\n $script_name ERROR: BAI $bai NOT GENERATED \n" >&2
	# delete BAM since something went wrong and it might be corrupted
	rm -fv "$bam"
	exit 1
fi


#########################


# run flagstat

bash_cmd="$sambamba_bin flagstat $bam > $flagstat_txt"
echo "CMD: $bash_cmd"
eval "$bash_cmd"


#########################


# check that flagstat output generated

if [ ! -s "$flagstat_txt" ] ; then
	echo -e "\n $script_name ERROR: FLAGSTAT $flagstat_txt NOT GENERATED \n" >&2
	# delete BAM and BAI since something went wrong and they might be corrupted
	rm -fv "$bam"
	rm -fv "$bai"
	exit 1
fi


#########################


# generate alignment summary

reads_filtered=$(samtools view -c "$bam")
echo "READS FILTERED: $reads_filtered"

reads_mapped_pct=$(echo "(${reads_mapped}/${reads_input})*100" | bc -l | cut -c 1-4)
reads_mapped_pct="${reads_mapped_pct}%"

reads_chrM_pct=$(echo "(${reads_chrM}/${reads_input})*100" | bc -l | cut -c 1-4)
reads_chrM_pct="${reads_chrM_pct}%"

reads_filtered_pct=$(echo "(${reads_filtered}/${reads_input})*100" | bc -l | cut -c 1-4)
reads_filtered_pct="${reads_filtered_pct}%"

# header for summary file
echo "#SAMPLE,TOTAL READS,MAPPED READS %,CHR M READS %,MAPPED READS MQ30 %" > "$summary_csv"

# summarize log file
echo "${sample},${reads_input},${reads_mapped_pct},${reads_chrM_pct},${reads_filtered_pct}" >> "$summary_csv"

sleep 5

# combine all sample summaries
cat ${summary_dir}/*.${segment_name}.csv | LC_ALL=C sort -t ',' -k1,1 | uniq > "${proj_dir}/summary.${segment_name}.csv"


#########################


# add sample and BAM to sample sheet
echo "${sample},${bam}" >> "$samples_csv"

sleep 5


#########################



# end
