#!/bin/bash


# Bismark methylation extractor


# script filename
script_name=$(basename "${BASH_SOURCE[0]}")
segment_name=${script_name/%.sh/}
echo -e "\n ========== $script_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 5 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name [project dir] [sample] [threads] [Bismark BAM] [analysis type] \n" >&2
	exit 1
fi

# arguments
proj_dir=$1
sample=$2
threads=$3
bismark_bam=$4
analysis_type=$5


#########################


# settings and files

segment_name="${segment_name}-${analysis_type}"

summary_dir="${proj_dir}/summary"
mkdir -p "$summary_dir"
summary_csv="${summary_dir}/${sample}.${segment_name}.csv"

bismark_meth_dir="${proj_dir}/meth-Bismark-${analysis_type}"
bismark_report_dir="${proj_dir}/Bismark-report"

bismark_report="${bismark_meth_dir}/${sample}_splitting_report.txt"
bismark_report_short="${bismark_meth_dir}/${sample}.report.txt"

bismark_mbias_report="${bismark_meth_dir}/${sample}.M-bias.txt"

bismark_bedgraph="${bismark_meth_dir}/${sample}.bedGraph"
bismark_bedgraph_gz="${bismark_bedgraph}.gz"

bismark_cpg_report_gz="${bismark_meth_dir}/${sample}.CpG_report.txt.gz"

# BISMARK_METHYLKIT_REPORT="${bismark_meth_dir}/${sample}.methylkit.txt"


#########################


# exit if output exits already

if [ -s "$bismark_report_short" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 1
fi


#########################


# check inputs and references

if [ ! -d "$proj_dir" ] ; then
	echo -e "\n $script_name ERROR: DIR $proj_dir DOES NOT EXIST \n" >&2
	exit 1
fi

if [ ! -s "$bismark_bam" ] ; then
	echo -e "\n $script_name ERROR: bismark_bam $bismark_bam DOES NOT EXIST \n" >&2
	exit 1
fi

code_dir=$(dirname "$(dirname "${BASH_SOURCE[0]}")")
ref_bismark=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" REF-BISMARK);

if [ ! -d "$ref_bismark" ] || [ ! "$ref_bismark" ] ; then
	echo -e "\n $script_name ERROR: REF $ref_bismark DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# adjust paramaters based on run type

# --ignore_r2 <int>
# The first couple of bases in Read 2 of BS-Seq experiments show a severe bias towards
# non-methylation as a result of end-repairing sonicated fragments with unmethylated cytosines,
# it is recommended that the first couple of bp of Read 2 are removed before starting downstream analysis.


if [ "$analysis_type" == "pe" ] ; then
	bismark_flags="--paired-end"
fi

# the first bases in R2 of BS-Seq experiments show a bias towards non-methylation as a result of end-repairing sonicated fragments with unmethylated cytosines
if [ "$analysis_type" == "pe-ignore-r2-3" ] ; then
	bismark_flags="--paired-end --ignore_r2 3"
fi

if [ "$analysis_type" == "se" ] ; then
	bismark_flags="--single-end"
fi

if [ "$analysis_type" == "se-ignore-r1-3" ] ; then
	bismark_flags="--single-end --ignore 3"
fi


#########################


# bismark_methylation_extractor

# load module (loads bowtie2/2.2.6 and samtools/1.3)
module load bismark/0.16.3

# navigate to output dir
mkdir -p $bismark_meth_dir
cd $bismark_meth_dir

# number of parallel instances of Bismark to run
multicore=$(( threads / 3 ))

echo " * BISMARK: $(readlink -f $(which bismark_methylation_extractor)) "
echo " * OUT DIR: $bismark_meth_dir "
echo " * BAM: $bismark_bam "
echo " * REPORT ORIGINAL: $bismark_report "
echo " * BEDGRAPH: $bismark_bedgraph_gz "
echo " * REPORT: $bismark_report_short "
echo " * CYTOSINE REPORT: $bismark_cpg_report_gz "

# --comprehensive \

# Genome-wide cytosine methylation report specific options:
# --cytosine_report - produces a genome-wide methylation report for all cytosines in the genome (1-based)
# --genome_folder - genome folder you wish to use to extract sequences from

bash_cmd="
bismark_methylation_extractor \
--gzip \
--multicore $multicore \
--buffer_size 4G \
--genome_folder $ref_bismark \
--cytosine_report \
--bedGraph \
$bismark_flags \
$bismark_bam
"
echo "CMD: $bash_cmd"
eval "$bash_cmd"


#########################


# check that output generated

if [ ! -s $bismark_report ] || [ ! $bismark_report ]
then
	printf "\n\n ERROR! REPORT $bismark_report NOT GENERATED \n\n"
	exit 1
fi


#########################


# put reports in report directory for bismark2report
# bismark_bams for bismark2summary (may not be necessary)

# adjust for SE/PE
if [ "$analysis_type" == "pe" ] ; then
	bismark_report_splittng="${bismark_report_dir}/${sample}_bismark_bt2_pe.deduplicated_splitting_report.txt"
	bismark_report_mbias="${bismark_report_dir}/${sample}_bismark_bt2_pe.deduplicated.M-bias.txt"
fi

# not tested
if [ "$analysis_type" == "se" ] ; then
	bismark_report_splittng="${bismark_report_dir}/${sample}_bismark_bt2_splitting_report.txt"
	bismark_report_mbias="${bismark_report_dir}/${sample}_bismark_bt2.M-bias.txt"
fi

if [ "$analysis_type" == "pe" ] || [ "$analysis_type" == "se" ]
then

	CMD="cp -v $bismark_report $bismark_report_splittng"
	printf "\n\n $CMD \n\n"
	$CMD

	CMD="cp -v $bismark_mbias_report $bismark_report_mbias"
	printf "\n\n $CMD \n\n"
	$CMD

fi


#########################


CMD="mv -v $bismark_report $bismark_report_short"
printf "\n\n $CMD \n\n"
$CMD


#########################


# convert cytosine report to methylkit-compatible format and remove uncovered bases

# cytosine report:
# <chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>
# new format: chr, start, strand, methylationRatio, coverage
# zcat $bismark_cpg_report_gz | awk '{ OFS="\t"; if( $4+0 > 0 || $5+0 > 0 ) print $1,$2,$3,$4/($4+$5),$4+$5; }' | sort -k1,1 -k2,2n > $BISMARK_METHYLKIT_REPORT


#########################


# convert bedgraph to bigwig

bigwig_dir="${proj_dir}/BIGWIG-Bismark"
bigwig="${bigwig_dir}/${sample}.bw"

mkdir -p "$bigwig_dir"

ref_chromsizes=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" REF-CHROMSIZES);

# bedgraph is gzipped by default
gunzip -cv $bismark_bedgraph_gz > $bismark_bedgraph

bash_cmd="bedGraphToBigWig $bismark_bedgraph $ref_chromsizes $bigwig"
echo "CMD: $bash_cmd"
eval "$bash_cmd"

# delete bedgraph
rm -fv $bismark_bedgraph


#########################


# generate summary

# header
echo "#SAMPLE,READS,TOTAL Cs,METH_CpG_CONTEXT,METH_CHG_CONTEXT,METH_CHH_CONTEXT" > $summary_csv

# print the relevant numbers
paste -d ',' \
<(echo "$sample") \
<(cat $bismark_report_short | grep "lines in total" 					| head -1 | cut -d " " -f 2) \
<(cat $bismark_report_short | grep "Total number of C's analysed" 		| head -1 | cut -f 2) \
<(cat $bismark_report_short | grep "C methylated in CpG context" 		| head -1 | cut -f 2) \
<(cat $bismark_report_short | grep "C methylated in CHG context" 		| head -1 | cut -f 2) \
<(cat $bismark_report_short | grep "C methylated in CHH context" 		| head -1 | cut -f 2) \
>> $summary_csv

sleep 30

# combine all sample summaries
cat ${summary_dir}/*.${segment_name}.csv | LC_ALL=C sort -t ',' -k1,1 | uniq \
> "${proj_dir}/summary.${segment_name}.csv"


#########################


# combine summary from each step

combined_summary_csv="${proj_dir}/summary-combined.bs-bismark-${analysis_type}.csv"

# combine all summaries
bash_cmd="
bash ${code_dir}/scripts/join-many.sh , X \
${proj_dir}/summary.fastq-clean.csv \
${proj_dir}/summary.fastq-trim-trimgalore.csv \
${proj_dir}/summary.fastq-trim-trimmomatic.csv \
${proj_dir}/summary.align-bismark.csv \
${proj_dir}/summary.align-dedup-bismark.csv \
${proj_dir}/summary.${segment_name}.csv \
> $combined_summary_csv
"
(eval $bash_cmd)


#########################


# combine charts into a single png

combined_png_2X="${proj_dir}/summary.${segment_name}.mbias.2x.png"
combined_png_4X="${proj_dir}/summary.${segment_name}.mbias.4x.png"

rm -f $combined_png_3X
rm -f $combined_png_4X

# -geometry +20+20 = 20px x and y padding
# -tile 4x = 4 images wide
# montage -geometry +20+20 -tile 4x ${bismark_meth_dir}/*M-bias*.png $COMBINED_PNG
montage -geometry +20+20 -tile 2x -label %t -pointsize 18 -font Nimbus-Sans-Regular ${bismark_meth_dir}/*M-bias*.png $combined_png_2X
montage -geometry +20+20 -tile 4x -label %t -pointsize 18 -font Nimbus-Sans-Regular ${bismark_meth_dir}/*M-bias*.png $combined_png_4X


#########################


# bismark html report

cd "$bismark_report_dir"

# sample report
bismark2report

# summary report
bismark2summary --basename '_summary'


#########################



# end
