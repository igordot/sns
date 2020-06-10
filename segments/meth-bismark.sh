#!/bin/bash


# Bismark methylation extractor


# script filename
script_path="${BASH_SOURCE[0]}"
script_name=$(basename "$script_path")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 5 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir sample_name num_threads BAM analysis_type \n" >&2
	if [ $# -gt 0 ] ; then echo -e "\n ARGS: $* \n" >&2 ; fi
	exit 1
fi

# arguments
proj_dir=$(readlink -f "$1")
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

bismark_cov="${bismark_meth_dir}/${sample}.bismark.cov"
bismark_cov_gz="${bismark_cov}.gz"

bismark_cpg_report_gz="${bismark_meth_dir}/${sample}.CpG_report.txt.gz"

# unload all loaded modulefiles
module purge
module add default-environment


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

code_dir=$(dirname $(dirname "$script_path"))

ref_bismark=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" REF-BISMARK);

if [ ! -d "$ref_bismark" ] || [ ! "$ref_bismark" ] ; then
	echo -e "\n $script_name ERROR: REF $ref_bismark DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# adjust paramaters based on run type

if [ "$analysis_type" == "pe" ] ; then
	bismark_flags="--paired-end"
fi

# --ignore_r2 <int>
# The first couple of bases in Read 2 of BS-Seq experiments show a severe bias towards
# non-methylation as a result of end-repairing sonicated fragments with unmethylated cytosines,
# it is recommended that the first couple of bp of Read 2 are removed before starting downstream analysis.
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

# bismark/0.22.1 loads bowtie2 and samtools (no version specified)
module add bismark/0.22.1

# navigate to the output dir
mkdir -p "$bismark_meth_dir"
cd "$bismark_meth_dir"

# number of parallel instances of Bismark to run
multicore=$(( threads / 3 ))

echo
echo " * Bismark: $(readlink -f $(which bismark)) "
echo " * Bismark version: $(bismark --version | grep -m 1 'Version' | tr -s '[:blank:]') "
echo " * output dir: $bismark_meth_dir "
echo " * BAM: $bismark_bam "
echo " * report original: $bismark_report "
echo " * BedGraph: $bismark_bedgraph_gz "
echo " * report: $bismark_report_short "
echo " * cytosine report: $bismark_cpg_report_gz "
echo

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
echo -e "\n CMD: $bash_cmd \n"
eval "$bash_cmd"


#########################


# check that output generated

if [ ! -s $bismark_report ] ; then
	printf "\n\n ERROR! REPORT $bismark_report NOT GENERATED \n\n"
	exit 1
fi

if [ ! -s $bismark_bedgraph_gz ] ; then
	printf "\n\n ERROR! BEDGRAPH $bismark_bedgraph_gz NOT GENERATED \n\n"
	exit 1
fi

if [ ! -s $bismark_cov_gz ] ; then
	printf "\n\n ERROR! COV $bismark_cov_gz NOT GENERATED \n\n"
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


# convert bedgraph to bigwig

module add ucscutils/374

bigwig_dir="${proj_dir}/BIGWIG-Bismark"
bigwig="${bigwig_dir}/${sample}.bw"

mkdir -p "$bigwig_dir"

ref_chromsizes=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" REF-CHROMSIZES);

echo " * bedGraphToBigWig: $(readlink -f $(which bedGraphToBigWig)) "
echo " * chrom sizes: $ref_chromsizes "
echo " * BedGraph gz: $bismark_bedgraph_gz "
echo " * cov gz: $bismark_cov_gz "
echo " * BedGraph: $bismark_bedgraph "
echo " * bigWig: $bigwig "

# bedgraph is gzipped by default
# gunzip -cv $bismark_bedgraph_gz | grep -v "track" | LC_ALL=C sort -k1,1 -k2,2n > $bismark_bedgraph

# bismark2bedGraph also writes out a coverage file (using 1-based genomic genomic coordinates):
# <chromosome> <start pos> <end pos> <methylation percentage> <count methylated> <count unmethylated>

# using coverage file instead of bedGraph to filter out low quality CpGs
gunzip -c $bismark_cov_gz \
| grep -v 'track' \
| awk -F $'\t' 'BEGIN {OFS=FS} ($5+$6 >= 5) {print $1,$2-1,$3,$4}' \
| LC_ALL=C sort -k1,1 -k2,2n \
> $bismark_bedgraph

bash_cmd="bedGraphToBigWig $bismark_bedgraph $ref_chromsizes $bigwig"
echo -e "\n CMD: $bash_cmd \n"
eval "$bash_cmd"

# delete bedgraph
rm -v $bismark_bedgraph


#########################


# generate summary

# header
header_report="READS,TOTAL Cs,METH_CpG_CONTEXT,METH_CHG_CONTEXT,METH_CHH_CONTEXT"
header_cpg_report="CpGs 1X,CpGs 10X,CpGs 100X"
echo "#SAMPLE,${header_report},${header_cpg_report}" > $summary_csv

# cytosine report:
# <chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>

# print the relevant numbers
paste -d ',' \
<(echo "$sample") \
<(cat "$bismark_report_short" | grep "lines in total"               | head -1 | cut -d " " -f 2) \
<(cat "$bismark_report_short" | grep "Total number of C's analysed" | head -1 | cut -f 2) \
<(cat "$bismark_report_short" | grep "C methylated in CpG context"  | head -1 | cut -f 2) \
<(cat "$bismark_report_short" | grep "C methylated in CHG context"  | head -1 | cut -f 2) \
<(cat "$bismark_report_short" | grep "C methylated in CHH context"  | head -1 | cut -f 2) \
<(gunzip -c "$bismark_cpg_report_gz" | awk -F $'\t' '$4 + $5 >= 1'   | wc -l) \
<(gunzip -c "$bismark_cpg_report_gz" | awk -F $'\t' '$4 + $5 >= 10'  | wc -l) \
<(gunzip -c "$bismark_cpg_report_gz" | awk -F $'\t' '$4 + $5 >= 100' | wc -l) \
>> $summary_csv

sleep 5

# combine all sample summaries
cat ${summary_dir}/*.${segment_name}.csv | LC_ALL=C sort -t ',' -k1,1 | uniq \
> "${proj_dir}/summary.${segment_name}.csv"

sleep 5


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
${proj_dir}/summary.bam-dedup-bismark.csv \
${proj_dir}/summary.${segment_name}.csv \
> $combined_summary_csv
"
(eval $bash_cmd)


#########################


# combine charts into a single png

combined_png_2X="${proj_dir}/summary.${segment_name}.mbias.2x.png"
combined_png_4X="${proj_dir}/summary.${segment_name}.mbias.4x.png"

rm "$combined_png_3X"
rm "$combined_png_4X"

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
