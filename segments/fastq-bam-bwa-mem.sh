#!/bin/bash


# run BWA-MEM


# script filename
script_name=$(basename "${BASH_SOURCE[0]}")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ $# -lt 4 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir sample_name num_threads FASTQ_R1 [FASTQ_R2] \n" >&2
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

code_dir=$(dirname "$(dirname "${BASH_SOURCE[0]}")")
ref_bwa=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" REF-BWA);

if [ ! -s "$ref_bwa" ] || [ ! -n "$ref_bwa" ] ; then
	echo -e "\n $script_name ERROR: REF $ref_bwa DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# settings and files

summary_dir="${proj_dir}/summary"
mkdir -p "$summary_dir"
summary_csv="${summary_dir}/${sample}.${segment_name}.csv"

samples_csv="${proj_dir}/samples.${segment_name}.csv"

bwa_bam_dir="${proj_dir}/BAM-BWA"
mkdir -p "$bwa_bam_dir"
bam="${bwa_bam_dir}/${sample}.bam"

bwa_logs_dir="${proj_dir}/logs-${segment_name}"
mkdir -p "$bwa_logs_dir"
bwa_flagstat="${bwa_logs_dir}/${sample}.flagstat.txt"


#########################


# exit if output exits already

if [ -s "$bam" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 1
fi


#########################


# BWA

module load bwa/0.7.13

sambamba_bin="/ifs/home/id460/bin/sambamba"

echo " * bwa: $(readlink -f $(which bwa)) "
echo " * bwa version: $(bwa 2>&1 | grep -m 1 'Version') "
echo " * sambamba: $(readlink -f $(which $sambamba_bin)) "
echo " * sambamba version: $($sambamba_bin 2>&1 | head -1) "
echo " * BWA REF: $ref_bwa "
echo " * FASTQ R1: $fastq_R1 "
echo " * FASTQ R2: $fastq_R2 "
echo " * BAM: $bam "

# step 1: align with BWA-MEM
# step 2: convert SAM to BAM and remove low quality reads
# step 3: sort BAM

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
--nthreads=${threads} \
--filter='mapping_quality>=10' \
--format=bam \
--compression-level=0 \
/dev/stdin \
| \
$sambamba_bin sort \
--nthreads=${threads} \
--memory-limit=16GB \
--out=${bam} \
/dev/stdin
"
echo "CMD: $bash_cmd"
eval "$bash_cmd"

sleep 30


#########################


# check that output generated

# check that bam file size above 10 kb
du -s $bam
bam_size=$(du -s --apparent-size $bam | cut -f 1)
if [ $bam_size -lt 10 ] ; then
	echo -e "\n $script_name ERROR: BAM $bam TOO SMALL \n" >&2
	exit 1
fi


#########################


# run flagstat

bash_cmd="$sambamba_bin flagstat $bam > $bwa_flagstat"
echo "CMD: $bash_cmd"
eval "$bash_cmd"


#########################


# generate alignment summary

# get number of input reads (not pairs)
fastq_lines=$(zcat "$fastq_R1" | wc -l)
if [ -n "$fastq_R2" ] ; then
	reads_input=$(echo "${fastq_lines}/2" | bc)
else
	reads_input=$(echo "${fastq_lines}/4" | bc)
fi

reads_mapped=$(cat "$bwa_flagstat" | grep -m 1 "mapped (" | cut -d ' ' -f 1)
reads_chimeric=$(cat "$bwa_flagstat" | grep -m 1 "mate mapped to a different chr" | cut -d ' ' -f 1)

reads_mapped_pct=$(echo "(${reads_mapped}/${reads_input})*100" | bc -l | cut -c 1-4)
reads_mapped_pct="${reads_mapped_pct}%"

reads_chimeric_pct=$(echo "(${reads_chimeric}/${reads_mapped})*100" | bc -l | cut -c 1-4)
reads_chimeric_pct="${reads_chimeric_pct}%"

# header for summary file
echo "#SAMPLE,INPUT READS,MAPPED READS (MQ10),MAPPED %,CHIMERIC %" > "$summary_csv"

# summarize log file
echo "${sample},${reads_input},${reads_mapped},${reads_mapped_pct},${reads_chimeric_pct}" >> "$summary_csv"

sleep 30

# combine all sample summaries
cat ${summary_dir}/*.${segment_name}.csv | LC_ALL=C sort -t ',' -k1,1 | uniq > "${proj_dir}/summary.${segment_name}.csv"


#########################


# add sample and BAM to sample sheet
echo "${sample},${bam}" >> "$samples_csv"

sleep 30

# sort and remove duplicates in place in sample sheet
LC_ALL=C sort -t ',' -k1,1 -u -o "$samples_csv" "$samples_csv"

sleep 30


#########################



# end



exit






#########################


# output
DD_DIR=${BASE_DIR}/BAM-DD
BAM_DD=${DD_DIR}/${ID}.dd.bam

mkdir -p $DD_DIR

print_msg "REMOVE DUPS - $ID"

CMD="bash ${HOME}/pipeline/wrappers/picard-mark-dups.sh $BASE_DIR $ID $BAM_BWA $BAM_DD false"
printf "\n\n $CMD \n\n"
$CMD

BAM=$BAM_DD


#########################


# check that BAM generated

if [ ! -s $BAM ]
then
	printf "\n\n ERROR! $BAM DOES NOT EXIST \n\n"
	exit 1
fi


#########################


print_msg "SAMTOOLS FLAGSTAT - $ID"

bash ${HOME}/pipeline/wrappers/samtools-flagstat.sh $BASE_DIR $ID $BAM


#########################


# Picard CollectInsertSizeMetrics

print_msg "INSERT SIZE METRICS - $ID"

# qsub-if-needed.sh <base dir> <expected out file> <analysis type> <cmd parameters>
bash ${HOME}/pipeline/helpers/qsub-if-needed.sh $BASE_DIR ${ID}.txt picard-insert-size $GENOME $BASE_DIR $ID $BAM $REF_FASTA


#########################


# try to find the bed file again in case it was placed in the directory after the run was started

if [ ! -s $BED ] || [ ! $BED ]
then
	BED=$(find $BASE_DIR -type f -name "*.bed")
fi

if [ ! -s $BED ] || [ ! $BED ]
then
	BED=0
fi


#########################


COV_DD_DIR=${BASE_DIR}/COV-DD
mkdir -p $COV_DD_DIR

if [ -f ${COV_DD_DIR}/${ID}.cov.ref.sample_summary ]
then

	print_msg "SKIP COVERAGE - $ID"

else

	print_msg "GATK COVERAGE - $ID"

	bash /ifs/home/id460/pipeline/wrappers/gatk-cov-all.sh $THREADS $REF_FASTA $BED $ID $BAM $COV_DD_DIR 

fi


#########################


# GATK requires several R libraries (using the default R 3.0.2 probably)
# install.packages("ggplot2")
# install.packages("gsalib")
# install.packages("reshape")

if [ $GENOME == "hg19" ]
then

	print_msg "GATK REALIGNMENT/RECALIBRATION (hg19)"

	bash /ifs/home/id460/pipeline/gatk-ra-rc-hg19.sh $ID $BAM $BASE_DIR $THREADS

elif [ $GENOME == "mm10" ]
then
	
	print_msg "GATK REALIGNMENT/RECALIBRATION (mm10)"

	bash /ifs/home/id460/pipeline/gatk-ra-rc-mm10.sh $ID $BAM $BASE_DIR $THREADS

elif [ $GENOME == "dm3" ]
then
	
	print_msg "GATK REALIGNMENT/RECALIBRATION (dm3)"

	bash /ifs/home/id460/pipeline/gatk-ra-rc-dm3.sh $ID $BAM $BASE_DIR $THREADS

	# realignment only (if chrM only)
	# bash /ifs/home/id460/pipeline/gatk-ra.sh $THREADS $GENOME $ID $BAM_DD $BASE_DIR

else

	print_msg "GENERIC GATK REALIGNMENT (NO RECALIBRATION, SINCE NO SNP INFO)"

	bash /ifs/home/id460/pipeline/gatk-ra.sh $THREADS $GENOME $ID $BAM $BASE_DIR

fi


#########################


# combine all summaries

# in case files do not exist (if re-running analysis, for example)
touch -a summary.trim.txt
touch -a summary.flagstat.txt

# (should change from paste to join)
# module load coreutils/8.24
paste \
${BASE_DIR}/summary.trim.txt \
${BASE_DIR}/summary.flagstat.txt \
${BASE_DIR}/COV-DD/__cov.txt \
${BASE_DIR}/COV-DD-RA-RC/_cov_samples.txt \
> ${BASE_DIR}/summary-combined.txt


#########################

date


# end


