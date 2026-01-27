#!/bin/bash


# run BWA-MEM


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

bwa_bam_dir="${proj_dir}/BAM-BWA"
mkdir -p "$bwa_bam_dir"
bam_temp="${bwa_bam_dir}/${sample}.tmp.bam"
bam_final="${bwa_bam_dir}/${sample}.bam"

bwa_logs_dir="${proj_dir}/logs-${segment_name}"
mkdir -p "$bwa_logs_dir"
bwa_flagstat="${bwa_logs_dir}/${sample}.flagstat.txt"

# unload all loaded modulefiles
module purge
module add default-environment


#########################


# check for output

# skip if final BAM and BAI exist
if [ -s "$bam_final" ] && [ -s "${bam_final}.bai" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	echo -e "\n $script_name ADD $sample TO $samples_csv \n" >&2
	echo "${sample},${bam_final}" >> "$samples_csv"
	exit 0
fi

# delete temp BAM (likely incomplete since the final BAM was not generated)
if [ -s "$bam_temp" ] ; then
	echo -e "\n $script_name WARNING: POTENTIALLY CORRUPT BAM $bam_temp EXISTS \n" >&2
	rm -fv "$bam_temp"
	rm -fv "${bam_temp}.bai"
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

ref_bwa=$(bash "${code_dir}/scripts/get-set-setting.sh" "${proj_dir}/settings.txt" REF-BWA)

if [ ! -s "$ref_bwa" ] || [ ! -n "$ref_bwa" ] ; then
	echo -e "\n $script_name ERROR: REF $ref_bwa DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# BWA

module add bwa/0.7.17
sambamba_bin="/gpfs/share/apps/sambamba/1.0.1/sambamba"

# check that the binary can be run
if ! bwa 2>&1 | grep -q "Version"; then
	echo -e "\n $script_name ERROR: bwa cannot be executed at $(which bwa) \n" >&2
	exit 1
fi
if ! $sambamba_bin --version >/dev/null 2>&1; then
	echo -e "\n $script_name ERROR: sambamba cannot be executed at $sambamba_bin \n" >&2
	exit 1
fi

echo
echo " * BWA: $(readlink -f $(which bwa)) "
echo " * BWA version: $(bwa 2>&1 | grep -m 1 'Version') "
echo " * sambamba: $(readlink -f $sambamba_bin) "
echo " * sambamba version: $($sambamba_bin 2>&1 | grep -m 1 'sambamba') "
echo " * BWA index: $ref_bwa "
echo " * FASTQ R1: $fastq_R1 "
echo " * FASTQ R2: $fastq_R2 "
echo " * BAM temp: $bam_temp "
echo " * BAM final: $bam_final "
echo

# step 1: align with BWA-MEM
# step 2: convert SAM to BAM and remove low quality reads
# step 3: sort BAM

# -M flag marks shorter split hits as secondary (for Picard compatibility)

bash_cmd="
bwa mem \
-M -v 1 \
-t $threads \
-R '@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA' \
$ref_bwa \
$fastq_R1 $fastq_R2 \
| \
$sambamba_bin view \
--sam-input \
--filter='mapping_quality>=10' \
--format=bam \
--compression-level=0 \
/dev/stdin \
| \
$sambamba_bin sort \
--memory-limit=8GB \
--out=${bam_temp} \
/dev/stdin
"
echo "CMD: $bash_cmd"
eval "$bash_cmd"

sleep 5


#########################


# check that bwa/sambamba output generated

if [ ! -s "$bam_temp" ] ; then
	echo -e "\n $script_name ERROR: BAM $bam_temp NOT GENERATED \n" >&2
	exit 1
fi

if [ ! -s "${bam_temp}.bai"] ; then
	echo -e "\n $script_name ERROR: BAI ${bam_dd_temp}.bai NOT GENERATED \n" >&2
	# delete BAM since something went wrong and it might be corrupted
	rm -fv "$bam_temp"
	exit 1
fi


#########################


# run flagstat

bash_cmd="$sambamba_bin flagstat $bam_temp > $bwa_flagstat"
echo "CMD: $bash_cmd"
eval "$bash_cmd"

sleep 5


#########################


# check that flagstat output generated

if [ ! -s "$bwa_flagstat" ] ; then
	echo -e "\n $script_name ERROR: FLAGSTAT $bwa_flagstat NOT GENERATED \n" >&2
	# delete BAM and BAI since something went wrong and they might be corrupted
	rm -fv "$bam_temp"
	rm -fv "${bam_temp}.bai"
	exit 1
fi


#########################


# generate alignment summary

# get number of input reads (not pairs)
fastq_lines=$(zcat "$fastq_R1" | wc -l)
if [ -n "$fastq_R2" ] ; then
	reads_input=$(echo "${fastq_lines}/2" | bc)
else
	reads_input=$(echo "${fastq_lines}/4" | bc)
fi

# "mapped" and "total" lines include secondary alignments
reads_mapped_all=$(cat "$bwa_flagstat" | grep -m 1 "mapped (" | cut -d ' ' -f 1)
reads_mapped_sec=$(cat "$bwa_flagstat" | grep -m 1 "secondary" | cut -d ' ' -f 1)
reads_mapped=$(echo "${reads_mapped_all} - ${reads_mapped_sec}" | bc)

reads_chimeric=$(cat "$bwa_flagstat" | grep -m 1 "mate mapped to a different chr" | cut -d ' ' -f 1)

reads_mapped_pct=$(echo "(${reads_mapped} / ${reads_input}) * 100" | bc -l | cut -c 1-4)
reads_mapped_pct="${reads_mapped_pct}%"

reads_chimeric_pct=$(echo "(${reads_chimeric} / ${reads_mapped}) * 100" | bc -l | cut -c 1-4)
reads_chimeric_pct="${reads_chimeric_pct}%"

# header for summary file
echo "#SAMPLE,INPUT READS,MAPPED READS (MQ10),MAPPED %,CHIMERIC %" > "$summary_csv"

# summarize log file
echo "${sample},${reads_input},${reads_mapped},${reads_mapped_pct},${reads_chimeric_pct}" >> "$summary_csv"

sleep 5

# combine all sample summaries
cat ${summary_dir}/*.${segment_name}.csv | LC_ALL=C sort -t ',' -k1,1 | uniq > "${proj_dir}/summary.${segment_name}.csv"


#########################


# finalize

# move BAM and BAI from temp to final
mv -v "$bam_temp" "${bam_final}"
mv -v "${bam_temp}.bai" "${bam_final}.bai"

# add sample and BAM to sample sheet
echo "${sample},${bam_final}" >> "$samples_csv"

sleep 5


#########################



# end
