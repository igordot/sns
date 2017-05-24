#!/bin/bash


# call copy number variants with Control-FREEC (WES settings, hg19 only)


# script filename
script_name=$(basename "${BASH_SOURCE[0]}")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 5 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir tumor_sample_name tumor_bam normal_sample_name normal_bam \n" >&2
	exit 1
fi

# arguments
proj_dir=$1
sample_t=$2
bam_t=$3
sample_n=$4
bam_n=$5


#########################


# check that inputs exist

if [ ! -d "$proj_dir" ] ; then
	echo -e "\n $script_name ERROR: PROJ DIR $proj_dir DOES NOT EXIST \n" >&2
	exit 1
fi

if [ ! -s "$bam_t" ] ; then
	echo -e "\n $script_name ERROR: BAM $bam_t DOES NOT EXIST \n" >&2
	exit 1
fi

if [ ! -s "$bam_n" ] ; then
	echo -e "\n $script_name ERROR: BAM $bam_n DOES NOT EXIST \n" >&2
	exit 1
fi

code_dir=$(dirname "$(dirname "${BASH_SOURCE[0]}")")

genome_dir=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" GENOME-DIR)

if [ ! -d "$genome_dir" ] ; then
	echo -e "\n $script_name ERROR: GENOME DIR $genome_dir DOES NOT EXIST \n" >&2
	exit 1
fi

ref_fasta=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" REF-FASTA)

if [ ! -s "$ref_fasta" ] ; then
	echo -e "\n $script_name ERROR: FASTA $ref_fasta DOES NOT EXIST \n" >&2
	exit 1
fi

found_probes_bed=$(find $proj_dir -maxdepth 1 -type f -iname "*probes*.bed" | head -1)
probes_bed_original=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" EXP-PROBES-BED "$found_probes_bed")

if [ ! -s "$probes_bed_original" ] ; then
	echo -e "\n $script_name ERROR: BED $probes_bed_original DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# settings and files

sample="${sample_t}:${sample_n}"

cnv_freec_dir="${proj_dir}/CNV-FREEC"
mkdir -p "$cnv_freec_dir"

# need a separate directory for each sample since some auto-generated files can have same filenames
sample_freec_logs_dir="${proj_dir}/logs-${segment_name}/${sample_t}-${sample_n}"
mkdir -p "$sample_freec_logs_dir"

config_txt="${sample_freec_logs_dir}/config.txt"
probes_bed_fixed="${sample_freec_logs_dir}/probes.bed"
chrs_txt="${sample_freec_logs_dir}/chrs.txt"
chrs_len_txt="${sample_freec_logs_dir}/chrs.len.txt"

out_base_sample="${sample_freec_logs_dir}/$(basename $bam_t)"
out_base_control="${sample_freec_logs_dir}/$(basename $bam_n)"
fixed_base="${cnv_freec_dir}/${sample_t}-${sample_n}"

cpn_sample="${out_base_sample}_sample.cpn"
cpn_control="${out_base_control}_control.cpn"

cnvs_original="${out_base_sample}_CNVs"
ratio_original="${out_base_sample}_ratio.txt"

cnvs_fixed="${fixed_base}.CNVs.txt"
ratio_fixed="${fixed_base}.ratio.txt"

graph_fixed="${fixed_base}.png"


#########################


# skip to annotation if output exits already

# annot_cmd="bash ${code_dir}/segments/annot-annovar.sh $proj_dir $sample $vcf_fixed"
# annot_cmd="bash ${code_dir}/segments/annot-regions-annovar.sh $proj_dir $sample $cnvs_fixed"

# if [ -s "$vcf_fixed" ] ; then
# 	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
# 	echo -e "\n CMD: $annot_cmd \n"
# 	($annot_cmd)
# 	exit 1
# fi

if [ -s "$cpn_sample" ] || [ -s "$cnvs_original" ] || [ -s "$cnvs_fixed" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 1
fi


#########################


# genome-specific settings (only hg19 supported for now)

genome_build=$(basename "$genome_dir")
if [[ "$genome_build" == "hg19" ]] ; then
	true
else
	echo -e "\n $script_name ERROR: UNSUPPORTED GENOME \n" >&2
	exit 1
fi


#########################


# generate FREEC-compatible references

module load bedtools/2.26.0

# clean up probes BED
fix_probes_cmd="
cat $probes_bed_original | cut -f 1-3 | grep -v '^chrM' | LC_ALL=C sort -u -k1,1 -k2,2n \
> $probes_bed_fixed
"
echo -e "\n CMD: $fix_probes_cmd \n"
eval "$fix_probes_cmd"

# list of chromosomes present in probes BED
chrs_txt_cmd="cat $probes_bed_fixed | cut -f 1 | uniq > $chrs_txt"
echo -e "\n CMD: $chrs_txt_cmd \n"
eval "$chrs_txt_cmd"

# chr lengths (must match BED file chrs, so using genome.fa.fai caused problems due to chrM)
chrs_len_txt_cmd="
cat ${ref_fasta}.fai | cut -f 1-2 | grep -F -w -f $chrs_txt > $chrs_len_txt
"
echo -e "\n CMD: $chrs_len_txt_cmd \n"
eval "$chrs_len_txt_cmd"


#########################


# create config

config_contents="

# whole exome sequencing config
# based on: https://github.com/BoevaLab/FREEC/blob/master/data/config_exome.txt


[general]

# output directory
outputDir = .

# number of threads
maxThreads = 4

# path to sambamba (used only to read .BAM files)
sambamba = /ifs/home/id460/software/sambamba/sambamba_v0.6.6

# file with chromosome lengths
# files of type hg19.fa.fai are also accepted starting from v9.3
chrLenFile = $chrs_len_txt

# path to the directory with chromosomes fasta files
# necessary to calculate a GC-content profile if a control dataset and GC-content profile are not available
chrFiles = /ifs/home/id460/ref/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/

# information about mappable positions (GEM output)
gemMappabilityFile = /ifs/home/id460/ref/hg19/FREEC/out100m2_hg19.gem

# genome ploidy
# you can set different values and Control-FREEC will select the one that explains most observed CNAs
ploidy = 2

# sample sex
# sex=XY will not annotate one copy of chr X and Y as a losssex=XY
sex = XY

# explicit window size
# for whole exome sequencing: window=0
window = 0

# set to 1 or 2 to correct the Read Count (RC) for GC-content bias and low mappability
# Default (WGS): 0
# Default (WES): 1 (≥ v9.5) and 0 (< v9.5)
# forceGCcontentNormalization = 1

# degree of polynomial
# Default: 3&4 (GC-content based normalization, WGS) or 1 (control-read-count-based normalization, WES)
# degree = 1

# desired behavior in the ambiguous regions
# 4: make a separate fragment of this unknown region and do not assign any copy number to this region at all
breakPointType = 4

# positive value of threshold for segmentation of normalized profiles
# Default: 0.8 (for WGS)
breakPointThreshold = 1.2

# set TRUE for target resequencing data (e.g., exome-seq)
noisyData = TRUE

# set FALSE to avoid printing -1 to the _ratio.txt files (useful for exome or targeted sequencing data)
printNA = FALSE

# threshold on the minimal number of reads per window in the control sample
# recommended value >=50 for for exome data
readCountThreshold = 50

# minimal number of consecutive windows to call a CNA
# Default: 3 for WES and 1 for WGS
minCNAlength = 5

# additional output in BedGraph format for the UCSC genome browser
BedGraphOutput = TRUE


[sample]

# file with mapped reads (can be single end reads, mate-pairs or paired-end reads)
mateFile = $bam_t

# format of reads (in mateFile)
# SAM, BAM, pileup, bowtie, eland, arachne, psl (BLAT), BED, Eland
inputFormat = BAM

# format of reads (in mateFile)
# 0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)
mateOrientation = FR


[control]

mateFile = $bam_n

inputFormat = BAM

mateOrientation = FR


[target]

# file with capture regions in .BED format
captureRegions = $probes_bed_fixed

"

echo "$config_contents" > "$config_txt"

sleep 30


#########################


# Control-FREEC

cd "$sample_freec_logs_dir"

freec_dir="/ifs/home/id460/software/FREEC/FREEC-10.6"
freec_bin="${freec_dir}/src/freec"

echo " * FREEC: $(readlink -f $freec_bin) "
echo " * sample T : $sample_t "
echo " * BAM T : $bam_t "
echo " * sample N : $sample_n "
echo " * BAM N : $bam_n "
echo " * probes original: $probes_bed_original "
echo " * probes fixed: $probes_bed_fixed "
echo " * CNVs original: $cnvs_original "
echo " * CNVs fixed: $cnvs_fixed "
echo " * ratio original: $ratio_original "
echo " * ratio fixed: $ratio_fixed "
echo " * graph: $graph_fixed "

freec_cmd="$freec_bin -conf $config_txt"
echo -e "\n CMD: $freec_cmd \n"
$freec_cmd

sleep 30


#########################


# check that output generated

if [ ! -s "$cnvs_original" ] ; then
	echo -e "\n $script_name ERROR: $cnvs_original NOT GENERATED \n" >&2
	exit 1
fi

if [ ! -s "$ratio_original" ] ; then
	echo -e "\n $script_name ERROR: $ratio_original NOT GENERATED \n" >&2
	exit 1
fi


#########################


# clean up

# delete raw copy number profiles
rm -fv "$cpn_sample"
rm -fv "$cpn_control"


#########################


# post-processing

module load r/3.3.0

# required libraries: rtracklayer

# add p-values (Wilcoxon test and Kolmogorov-Smirnov test) and columns header to the predicted CNVs
freec_asses_sig_cmd="cat ${freec_dir}/scripts/assess_significance.R | R --slave --args $cnvs_original $ratio_original"
echo -e "\n CMD: $freec_asses_sig_cmd \n"
eval "$freec_asses_sig_cmd"

# add "chr" to CNV table chromosomes names
cat ${cnvs_original}.p.value.txt | sed 's/^\([0-9XY]\)/chr\1/' | LC_ALL=C sort -k1,1 -k2,2n | uniq > $cnvs_fixed

# visualize normalized copy number profile with predicted CNAs as well as BAF profile by running makeGraph.R
freec_makegraph_cmd="cat ${freec_dir}/scripts/makeGraph.R | R --slave --args 2 $ratio_original"
echo -e "\n CMD: $freec_makegraph_cmd \n"
eval "$freec_makegraph_cmd"


#########################


# fix some of the names

# mv -v "$cnvs_original" "$cnvs_fixed"
mv -v "$ratio_original" "$ratio_fixed"
mv -v "${ratio_original}.png" "$graph_fixed"

sleep 30


#########################


# check that output generated

if [ ! -s "$cnvs_fixed" ] ; then
	echo -e "\n $script_name ERROR: CNVs $cnvs_fixed NOT GENERATED \n" >&2
	exit 1
fi


#########################


# annotate

annot_cmd="bash ${code_dir}/segments/annot-regions-annovar.sh $proj_dir $sample $cnvs_fixed"
echo -e "\n CMD: $annot_cmd \n"
($annot_cmd)


#########################



# end
