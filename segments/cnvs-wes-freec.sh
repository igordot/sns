#!/bin/bash


# call copy number variants with Control-FREEC (WES settings)


# script filename
script_path="${BASH_SOURCE[0]}"
script_name=$(basename "$script_path")
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

code_dir=$(dirname $(dirname "$script_path"))

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

chrom_sizes=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" REF-CHROMSIZES)

if [ ! -s "$chrom_sizes" ] ; then
	echo -e "\n $script_name ERROR: CHROM SIZES $chrom_sizes DOES NOT EXIST \n" >&2
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

summary_dir="${proj_dir}/summary"
mkdir -p "$summary_dir"
summary_csv="${summary_dir}/${sample}.${segment_name}.csv"

cnv_freec_dir="${proj_dir}/CNV-FREEC"
mkdir -p "$cnv_freec_dir"

# need a separate directory for each sample since some auto-generated files can have same filenames
sample_freec_logs_dir="${proj_dir}/logs-${segment_name}/${sample_t}-${sample_n}"
mkdir -p "$sample_freec_logs_dir"

config_txt="${sample_freec_logs_dir}/config.txt"
probes_bed_fixed="${sample_freec_logs_dir}/probes.bed"
chrs_txt="${sample_freec_logs_dir}/chrs.txt"
chrs_len_txt="${sample_freec_logs_dir}/chrs.len.txt"
snps_bed_fixed="${sample_freec_logs_dir}/snps.bed"

out_base_sample="${sample_freec_logs_dir}/$(basename $bam_t)"
out_base_control="${sample_freec_logs_dir}/$(basename $bam_n)"
fixed_base="${cnv_freec_dir}/${sample_t}-${sample_n}"

cpn_sample="${out_base_sample}_sample.cpn"
cpn_control="${out_base_control}_control.cpn"

minipileup_sample="${out_base_sample}_minipileup.pileup"
minipileup_control="${out_base_control}_minipileup.pileup"

cnvs_original="${out_base_sample}_CNVs"
ratio_original="${out_base_sample}_ratio.txt"
bafs_original="${out_base_sample}_BAF.txt"

cnvs_fixed="${fixed_base}.CNVs.txt"
ratio_fixed="${fixed_base}.ratio.txt"
bafs_fixed="${fixed_base}.BAFs.txt"

graph_fixed="${fixed_base}.png"

# annotation command (next segment)
annot_cmd="bash ${code_dir}/segments/annot-regions-annovar.sh $proj_dir $sample $cnvs_fixed"

# unload all loaded modulefiles
module purge
module add default-environment


#########################


# check for output

# skip to annotation if final output exists already
if [ -s "$graph_fixed" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	echo -e "\n CMD: $annot_cmd \n"
	($annot_cmd)
	exit 0
fi

# delete original CNVs file (likely incomplete since the final graph was not generated)
if [ -s "$cnvs_original" ] ; then
	echo -e "\n $script_name WARNING: POTENTIALLY CORRUPT $cnvs_original EXISTS \n" >&2
	rm -fv "$cnvs_original"
	# delete all output
	rm -rf "${sample_freec_logs_dir}"
fi

# delete original cpn file (likely incomplete since the final graph was not generated)
if [ -s "$cpn_sample" ] ; then
	echo -e "\n $script_name WARNING: POTENTIALLY CORRUPT $cpn_sample EXISTS \n" >&2
	rm -fv "$cpn_sample"
	# delete all output
	rm -rf "${sample_freec_logs_dir}"
fi


#########################


# genome-specific settings

genome_build=$(basename "$genome_dir")
if [[ "$genome_build" == "hg19" ]] ; then
	chr_files_dir="/gpfs/data/igorlab/ref/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/"
	gem="/gpfs/data/igorlab/ref/hg19/FREEC/out100m2_hg19.gem"
	snps_vcf="/gpfs/data/igorlab/ref/hg19/dbSNP/common_all_20170710.snv.maf5.vcf"
	snps_bed="/gpfs/data/igorlab/ref/hg19/dbSNP/common_all_20170710.snv.maf5.bed"
elif [[ "$genome_build" == "hg38" ]] ; then
	chr_files_dir="/gpfs/data/igorlab/ref/hg38/chromosomes/"
	gem="/gpfs/data/igorlab/ref/hg38/genome.len100.mm2.mappability"
	snps_vcf="/gpfs/data/igorlab/ref/hg38/dbSNP/common_all_20170710.snv.maf5.vcf"
	snps_bed="/gpfs/data/igorlab/ref/hg38/dbSNP/common_all_20170710.snv.maf5.bed"
elif [[ "$genome_build" == "mm10" ]] ; then
	chr_files_dir="/gpfs/data/igorlab/ref/iGenomes/Mus_musculus/UCSC/mm10/Sequence/Chromosomes/"
	gem="/gpfs/data/igorlab/ref/mm10/FREEC/out100m4_mm10.gem"
	snps_vcf=""
	snps_bed=""
	echo -e "\n $script_name ERROR: UNSUPPORTED GENOME \n" >&2
	exit 1
elif [[ "$genome_build" == "canFam3" ]] ; then
	chr_files_dir="/gpfs/data/igorlab/ref/canFam3/chromosomes/"
	gem="/gpfs/data/igorlab/ref/canFam3/genome.len100.mm2.mappability"
	snps_vcf="/gpfs/data/igorlab/ref/canFam3/dbSNP/dbsnp.151.snp.validated.vcf"
	snps_bed="/gpfs/data/igorlab/ref/canFam3/dbSNP/dbsnp.151.snp.validated.bed"
else
	echo -e "\n $script_name ERROR: UNSUPPORTED GENOME \n" >&2
	exit 1
fi

if [ ! -s "$gem" ] ; then
	echo -e "\n $script_name ERROR: GEM $gem DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# generate FREEC-compatible references

# FREEC compiled with GCC 6.1.0 (load same GCC when running)
module add gcc/6.1.0
# bedtools to create .pileup files for WES data
module add bedtools/2.27.1
# samtools to create .pileup files (for BAF) (even with sambamba enabled)
module add samtools/1.3

# clean up probes BED
fix_probes_cmd="
cat $probes_bed_original | cut -f 1-3 | grep -v '^chrM' | grep -v '^chrUn' | LC_ALL=C sort -u -k1,1 -k2,2n \
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

# filter SNPs BED file to only include positions covered by probes +/- 500bp
fix_snps_bed_cmd="
bedtools slop -i $probes_bed_fixed -g $chrom_sizes -b 500 \
| bedtools merge -i stdin -d 10 \
| bedtools intersect -wa -a $snps_bed -b stdin \
> $snps_bed_fixed
"
echo -e "\n CMD: $fix_snps_bed_cmd \n"
eval "$fix_snps_bed_cmd"

sleep 5


#########################


# create config

# whole exome sequencing config
# based on: https://github.com/BoevaLab/FREEC/blob/master/data/config_exome.txt

config_contents="

[general]

# output directory
outputDir = .

# number of threads
maxThreads = 4

# file with chromosome lengths
# files of type hg19.fa.fai are also accepted starting from v9.3
chrLenFile = $chrs_len_txt

# path to the directory with chromosomes fasta files
# necessary to calculate a GC-content profile if a control dataset and GC-content profile are not available
chrFiles = $chr_files_dir

# information about mappable positions (GEM output)
gemMappabilityFile = $gem

# genome ploidy
# you can set different values and Control-FREEC will select the one that explains most observed CNAs
ploidy = 2

# sample sex
# sex=XY will not annotate one copy of chr X and Y as a losssex=XY
sex = XY

# explicit window size
# for whole exome sequencing: window=0
window = 0

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
readCountThreshold = 25

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


[BAF]

# file with known SNPs (.vcf or .vcf.gz files are also accepted in >v9.3)
SNPfile = $snps_vcf

# BED or VCF file with SNP positions to create a pileup file from the BAM file provided in mateFile
makePileup = $snps_bed_fixed

# one fasta file for the whole genome (need to be provided only when makePileup is used)
fastaFile = $ref_fasta

# minimal read coverage for a position to be considered in BAF analysis (default 0)
minimalCoveragePerPosition = 10

# basis for Phred quality (Default: 0; usually 33 or 64)
shiftInQuality = 33


[target]

# file with capture regions in .BED format
captureRegions = $probes_bed_fixed

"

echo "$config_contents" > "$config_txt"

sleep 5


#########################


# Control-FREEC

cd "$sample_freec_logs_dir"

freec_dir="/gpfs/data/igorlab/software/FREEC/FREEC-11.6"
freec_bin="${freec_dir}/src/freec"

echo
echo " * FREEC: $(readlink -f $freec_bin) "
echo " * sample T : $sample_t "
echo " * BAM T : $bam_t "
echo " * sample N : $sample_n "
echo " * BAM N : $bam_n "
echo " * probes original: $probes_bed_original "
echo " * probes fixed: $probes_bed_fixed "
echo " * ratio original: $ratio_original "
echo " * ratio fixed: $ratio_fixed "
echo " * CNVs original: $cnvs_original "
echo " * CNVs fixed: $cnvs_fixed "
echo " * BAFs original: $bafs_original "
echo " * BAFs fixed: $bafs_fixed "
echo " * graph: $graph_fixed "
echo

freec_cmd="$freec_bin -conf $config_txt"
echo -e "\n CMD: $freec_cmd \n"
($freec_cmd)

sleep 5


#########################


# check that output generated

if [ ! -s "$ratio_original" ] ; then
	echo -e "\n $script_name ERROR: $ratio_original NOT GENERATED \n" >&2
	exit 1
fi

if [ ! -s "$cnvs_original" ] ; then
	echo -e "\n $script_name ERROR: $cnvs_original NOT GENERATED \n" >&2
	exit 1
fi

if [ ! -s "$bafs_original" ] ; then
	echo -e "\n $script_name ERROR: $bafs_original NOT GENERATED \n" >&2
	exit 1
fi

if [ ! -s "$minipileup_sample" ] ; then
	echo -e "\n $script_name ERROR: $minipileup_sample NOT GENERATED \n" >&2
fi

if [ ! -s "$minipileup_control" ] ; then
	echo -e "\n $script_name ERROR: $minipileup_control NOT GENERATED \n" >&2
fi


#########################


# clean up

# delete raw copy number profiles
rm -fv "$cpn_sample"
rm -fv "$cpn_control"

# delete minipileups
rm -fv "$minipileup_sample"
rm -fv "$minipileup_control"


#########################


# post-processing

# clean up the environment before loading R module to avoid GCC conflicts
module purge
module add default-environment
module add r/3.6.1

echo
echo " * R: $(readlink -f $(which R)) "
echo " * R version: $(R --version | head -1) "
echo " * Rscript: $(readlink -f $(which Rscript)) "
echo " * Rscript version: $(Rscript --version 2>&1) "
echo

# test relevant R packages
Rscript --vanilla "${code_dir}/scripts/test-package.R" rtracklayer

# add p-values (Wilcoxon test and Kolmogorov-Smirnov test) and columns header to the predicted CNVs
freec_asses_sig_cmd="cat ${freec_dir}/scripts/assess_significance.R | R --slave --args $cnvs_original $ratio_original"
echo -e "\n CMD: $freec_asses_sig_cmd \n"
eval "$freec_asses_sig_cmd"

# add "chr" to CNV table chromosomes names
cat "${cnvs_original}.p.value.txt" | sed 's/^\([0-9XY]\)/chr\1/' | LC_ALL=C sort -k1,1 -k2,2n | uniq > "$cnvs_fixed"

# add "chr" to BAF table chromosomes names
cat "$bafs_original" | sed 's/^\([0-9XY]\)/chr\1/' | LC_ALL=C sort -k1,1 -k2,2n | uniq > "$bafs_fixed"

# visualize normalized copy number profile with predicted CNAs as well as BAF profile by running makeGraph.R
freec_makegraph_cmd="cat ${freec_dir}/scripts/makeGraph.R | R --slave --args 2 $ratio_original"
echo -e "\n CMD: $freec_makegraph_cmd \n"
eval "$freec_makegraph_cmd"


#########################


# fix some of the names

mv -v "$ratio_original" "$ratio_fixed"
mv -v "${ratio_original}.png" "$graph_fixed"


#########################


# check that output generated

if [ ! -s "$cnvs_fixed" ] ; then
	echo -e "\n $script_name ERROR: CNVs $cnvs_fixed NOT GENERATED \n" >&2
	exit 1
fi


#########################


# summary

# ratios and predicted copy number alterations for each window
num_bins=$(cat "$ratio_fixed" | grep -v 'MedianRatio' | wc -l)
echo "num bins: $num_bins"

# B-allele frequencies for each possibly heterozygous SNP position
num_bafs=$(cat "$bafs_fixed" | grep -v 'FittedA' | wc -l)
echo "num BAF SNPs: $num_bafs"

# predicted copy number alterations
num_cnas=$(cat "$cnvs_fixed" | grep -v 'uncertainty' | wc -l)
echo "num CNAs: $num_cnas"

# header for summary file
echo "#SAMPLE,bins,BAF SNPs,CNAs" > "$summary_csv"

# summarize log file
echo "${sample},${num_bins},${num_bafs},${num_cnas}" >> "$summary_csv"

sleep 5

# combine all sample summaries
cat ${summary_dir}/*.${segment_name}.csv | LC_ALL=C sort -t ',' -k1,1 | uniq > "${proj_dir}/summary.${segment_name}.csv"


#########################


# annotate

echo -e "\n CMD: $annot_cmd \n"
($annot_cmd)


#########################



# end
