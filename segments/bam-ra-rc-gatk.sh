#!/bin/bash


# GATK realignment and recalibration


# script filename
script_name=$(basename "${BASH_SOURCE[0]}")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 4 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir sample_name num_threads BAM \n" >&2
	exit 1
fi

# arguments
proj_dir=$1
sample=$2
threads=$3
bam=$4


#########################


# check that inputs exist

if [ ! -d "$proj_dir" ] ; then
	echo -e "\n $script_name ERROR: PROJ DIR $proj_dir DOES NOT EXIST \n" >&2
	exit 1
fi

if [ ! -s "$bam" ] ; then
	echo -e "\n $script_name ERROR: BAM $bam DOES NOT EXIST \n" >&2
	exit 1
fi

code_dir=$(dirname "$(dirname "${BASH_SOURCE[0]}")")

genome_dir=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" GENOME-DIR);

if [ ! -d "$genome_dir" ] ; then
	echo -e "\n $script_name ERROR: GENOME DIR $genome_dir DOES NOT EXIST \n" >&2
	exit 1
fi

ref_fasta=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" REF-FASTA);

if [ ! -s "$ref_fasta" ] ; then
	echo -e "\n $script_name ERROR: FASTA $ref_fasta DOES NOT EXIST \n" >&2
	exit 1
fi

ref_dict=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" REF-DICT);

if [ ! -s "$ref_dict" ] ; then
	echo -e "\n $script_name ERROR: DICT $ref_dict DOES NOT EXIST \n" >&2
	exit 1
fi

found_bed=$(find $proj_dir -maxdepth 1 -type f -name "*.bed" | head -1)
bed=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" EXP-BED $found_bed);

if [ ! -s "$bed" ] ; then
	echo -e "\n $script_name ERROR: BED $bed DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# settings and files

samples_csv="${proj_dir}/samples.${segment_name}.csv"

# account for both dedup (.dd.bam) non-dedup (.bam) input BAMs
bam_base=$(basename "$bam")
bam_base=${bam_base/%.bam/}

bam_ra_dir="${proj_dir}/BAM-GATK-RA"
mkdir -p "$bam_ra_dir"
bam_ra="${bam_ra_dir}/${bam_base}.ra.bam"
bai_ra="${bam_ra_dir}/${bam_base}.ra.bai"

bam_ra_rc_dir="${proj_dir}/BAM-GATK-RA-RC"
mkdir -p "$bam_ra_rc_dir"
bam_ra_rc="${bam_ra_rc_dir}/${bam_base}.ra.rc.bam"
bai_ra_rc="${bam_ra_rc_dir}/${bam_base}.ra.rc.bai"

gatk_logs_dir="${proj_dir}/logs-${segment_name}"
mkdir -p "$gatk_logs_dir"
gatk_ra_intervals="${gatk_logs_dir}/${sample}.intervals"
gatk_rc_table1="${gatk_logs_dir}/${sample}.table1.txt"
gatk_rc_table2="${gatk_logs_dir}/${sample}.table2.txt"
gatk_rc_csv="${gatk_logs_dir}/${sample}.csv"
gatk_rc_pdf="${gatk_logs_dir}/${sample}.pdf"


#########################


# exit if output exits already

if [ -s "$bam_ra_rc" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 1
fi

if [ -s "$bam_ra" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 1
fi


#########################


# GATK settings

module unload java
module load java/1.8
module load r/3.3.0

# command
gatk_jar="/ifs/home/id460/software/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar"
gatk_cmd="java -Xms16G -Xmx16G -jar ${gatk_jar}"

if [ ! -s "$gatk_jar" ] ; then
	echo -e "\n $script_name ERROR: GATK $gatk_jar DOES NOT EXIST \n" >&2
	exit 1
fi

# error log (DEBUG, INFO (default), WARN, ERROR, FATAL, OFF)
gatk_log_level_arg="--logging_level ERROR"

# known variants (may vary greatly for each genome)
if [[ "$genome_dir" == */hg19 ]] ; then
	gatk_indel_vcf_1="${genome_dir}/gatk-bundle/1000G_phase1.indels.hg19.vcf"
	gatk_indel_vcf_2="${genome_dir}/gatk-bundle/Mills_and_1000G_gold_standard.indels.hg19.vcf"
	gatk_snp_vcf="${genome_dir}/gatk-bundle/dbsnp_138.hg19.vcf"
	gatk_ra_known_arg="-known $gatk_indel_vcf_1 -known $gatk_indel_vcf_2"
	gatk_rc_known_arg="-knownSites $gatk_indel_vcf_1 -knownSites $gatk_indel_vcf_2 -knownSites $gatk_snp_vcf"
elif [[ "$genome_dir" == */mm10 ]] ; then
	gatk_indel_vcf="${genome_dir}/MGP/mgp.v5.indels.pass.chr.sort.vcf"
	gatk_snp_vcf="${genome_dir}/dbSNP/dbsnp.146.vcf"
	gatk_ra_known_arg="-known $gatk_indel_vcf"
	gatk_rc_known_arg="-knownSites $gatk_indel_vcf -knownSites $gatk_snp_vcf"
elif [[ "$genome_dir" == */dm3$ ]] ; then
	gatk_ra_known_arg="-known xxxxx"
else
	# should adjust
	echo -e "\n $script_name ERROR: genome $genome_dir not supported \n" >&2
	exit 1
fi


#########################


# realignment

echo " * GATK: $(readlink -f $gatk_jar) "
echo " * GATK version: $($gatk_cmd --version) "
echo " * BAM IN: $bam "
echo " * BAM RA: $bam_ra "

gatk_ra1_cmd="
$gatk_cmd -T RealignerTargetCreator -dt NONE $gatk_log_level_arg \
-nt $threads \
--reference_sequence $ref_fasta \
$gatk_ra_known_arg \
--intervals $bed --interval_padding 10 \
--input_file $bam \
--out $gatk_ra_intervals
"
echo "CMD: $gatk_ra1_cmd"
$gatk_ra1_cmd

gatk_ra2_cmd="
$gatk_cmd -T IndelRealigner -dt NONE $gatk_log_level_arg \
--reference_sequence $ref_fasta \
--maxReadsForRealignment 50000 \
$gatk_ra_known_arg \
-targetIntervals $gatk_ra_intervals \
--input_file $bam \
--out $bam_ra
"
echo "CMD: $gatk_ra2_cmd"
$gatk_ra2_cmd

sleep 30


#########################


# check that output generated

if [ ! -s "$bam_ra" ] ; then
	echo -e "\n $script_name ERROR: BAM $bam_ra NOT GENERATED \n" >&2
	exit 1
fi


#########################


# recalibration

# first pass recalibration table file
gatk_rc1_cmd="
$gatk_cmd -T BaseRecalibrator $gatk_log_level_arg \
-nct $threads \
-rf BadCigar \
--reference_sequence $ref_fasta \
$gatk_rc_known_arg \
--intervals $bed --interval_padding 10 \
--input_file $bam_ra \
--out $gatk_rc_table1
"
echo "CMD: $gatk_rc1_cmd"
$gatk_rc1_cmd

# second pass recalibration table file
gatk_rc2_cmd="
$gatk_cmd -T BaseRecalibrator $gatk_log_level_arg \
-nct $threads \
-rf BadCigar \
--reference_sequence $ref_fasta \
$gatk_rc_known_arg \
--intervals $bed --interval_padding 10 \
--input_file $bam_ra \
-BQSR $gatk_rc_table1 \
--out $gatk_rc_table2
"
echo "CMD: $gatk_rc2_cmd"
$gatk_rc2_cmd

# generate the plots report and also keep a copy of the csv (optional)
gatk_rc3_cmd="
$gatk_cmd -T AnalyzeCovariates $gatk_log_level_arg \
--reference_sequence $ref_fasta \
-before $gatk_rc_table1 \
-after $gatk_rc_table2 \
-csv $gatk_rc_csv \
-plots $gatk_rc_pdf
"
echo "CMD: $gatk_rc3_cmd"
$gatk_rc3_cmd

# generate recalibrated BAM
gatk_rc4_cmd="
$gatk_cmd -T PrintReads $gatk_log_level_arg \
-nct $threads \
-rf BadCigar \
--reference_sequence $ref_fasta \
-BQSR $gatk_rc_table1 \
--input_file $bam_ra \
--out $bam_ra_rc
"
echo "CMD: $gatk_rc4_cmd"
$gatk_rc4_cmd

sleep 30


#########################


# check that output generated

if [ ! -s "$bam_ra_rc" ] ; then
	echo -e "\n $script_name ERROR: BAM $bam_ra_rc NOT GENERATED \n" >&2
	exit 1
fi


#########################


# clean up

# move .bai to .bam.bai (some tools expect that)
mv -v "$bai_ra_rc" "${bam_ra_rc}.bai"

# delete files that are no longer needed
rm -fv "$bam_ra"
rm -fv "$bai_ra"
rm -fv "$gatk_ra_intervals"
rm -fv "$gatk_rc_table1"
rm -fv "$gatk_rc_table2"


#########################


# add sample and BAM to sample sheet
echo "${sample},${bam_ra_rc}" >> "$samples_csv"

sleep 30

# sort and remove duplicates in place in sample sheet
LC_ALL=C sort -t ',' -k1,1 -u -o "$samples_csv" "$samples_csv"

sleep 30


#########################



# end
