#!/bin/bash


# annotate regions with cytoband and overlapping genes using ANNOVAR
# expects tab-delimited input file with header row and 5 columns: chr, start, end, 0, 0


# script filename
script_path="${BASH_SOURCE[0]}"
script_name=$(basename "$script_path")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 3 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir sample_name regions_table \n" >&2
	exit 1
fi

# arguments
proj_dir=$1
sample=$2
regions_table=$3


#########################


# settings and files

input_type=$(basename "$(dirname "$regions_table")")

segment_name="${input_type}-annot"

# summary_dir="${proj_dir}/summary"
# mkdir -p "$summary_dir"
# summary_csv="${summary_dir}/${sample}.${segment_name}.csv"

annovar_dir="${proj_dir}/${segment_name}"
mkdir -p "$annovar_dir"

# clean up sample name for paired samples
sample_clean=${sample//:/-}

annovar_input="${annovar_dir}/${sample_clean}.avinput"
annovar_out_prefix="${annovar_dir}/${sample_clean}"
annovar_out_fixed="${annovar_out_prefix}.annot.txt"
annovar_combined="${annovar_out_prefix}.combined.txt"
regions_table_fixed="${annovar_out_prefix}.in.txt"

# unload all loaded modulefiles
module purge
module load local


#########################


# exit if output exits already

if [ -s "$annovar_combined" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 1
fi


#########################


# check inputs and references

if [ ! -d "$proj_dir" ] ; then
	echo -e "\n $script_name ERROR: DIR $proj_dir DOES NOT EXIST \n" >&2
	exit 1
fi

if [ ! -s "$regions_table" ] ; then
	echo -e "\n $script_name ERROR: REGIONS TABLE $regions_table DOES NOT EXIST \n" >&2
	exit 1
fi

code_dir=$(dirname $(dirname "$script_path"))

genome_dir=$(bash "${code_dir}/scripts/get-set-setting.sh" "${proj_dir}/settings.txt" GENOME-DIR);

if [ ! -d "$genome_dir" ] ; then
	echo -e "\n $script_name ERROR: GENOME DIR $genome_dir DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# ANNOVAR genome-specific settings

# ANNOVAR directory
annovar_path="/ifs/home/id460/software/annovar/annovar-170716"
annovar_db_path="/ifs/home/id460/ref/annovar"

genome_build=$(basename "$genome_dir")

# table_annovar output file (automatically named)
annovar_multianno="${annovar_out_prefix}.${genome_build}_multianno.txt"

# genome-specific settings (available annotations differ)
# annovar_buildver, annovar_protocol, annovar_operation - table_annovar parameters
# annovar_keep_cols - columns to keep for final fixed table
if [[ "$genome_build" == "hg19" ]] ; then
	annovar_buildver="hg19"
	annovar_protocol="cytoBand,refGene"
	annovar_operation="r,g"
	annovar_keep_cols="1,7-10"
elif [[ "$genome_build" == "mm10" ]] ; then
	annovar_buildver="mm10"
	annovar_protocol="cytoBand,refGene"
	annovar_operation="r,g"
	annovar_keep_cols="1,7-10"
elif [[ "$genome_build" == "canFam3" ]] ; then
	annovar_buildver="canFam3"
	annovar_protocol="refGene"
	annovar_operation="g"
	annovar_keep_cols="1,8-10"
else
	annovar_buildver=""
	echo -e "\n $script_name ERROR: UNKNOWN GENOME $genome_build \n" >&2
	exit 1
fi


#########################


# ANNOVAR convert2annovar - convert VCF to ANNOVAR input format

echo " * convert2annovar path: $(readlink -f ${annovar_path}/convert2annovar.pl) "
echo " * regions table : $regions_table "
echo " * ANNOVAR out dir: $annovar_dir "
echo " * convert2annovar out : $annovar_input "

convert_cmd="
cat $regions_table \
| grep -v 'start.*end' \
| cut -f 1-3 \
| sort -k1,1 -k2,2n \
| uniq \
| awk -F $'\t' 'BEGIN {OFS=FS} {print \$0,\"0\",\"0\"}' \
> "$annovar_input"
"
echo -e "\n CMD: $convert_cmd \n"
eval "$convert_cmd"

sleep 30


#########################


# check that convert2annovar completed

if [ ! -s "$annovar_input" ] ; then
	echo -e "\n $script_name ERROR: $annovar_input IS EMPTY \n" >&2
	exit 1
fi


#########################


# ANNOVAR table_annovar - run a pipeline on a list of variants and summarize their functional effects

echo " * table_annovar path: $(readlink -f ${annovar_path}/table_annovar.pl) "
echo " * ANNOVAR out dir: $annovar_dir "
echo " * table_annovar out prefix : $annovar_out_prefix "
echo " * table_annovar out : $annovar_multianno "

# annotate with annovar (outputs $annovar_multianno)
table_cmd="
perl ${annovar_path}/table_annovar.pl $annovar_input ${annovar_db_path}/${annovar_buildver}/ \
--outfile $annovar_out_prefix \
--buildver $annovar_buildver \
--protocol $annovar_protocol \
--operation $annovar_operation \
--nastring . \
--remove
"
echo -e "\n CMD: $table_cmd \n"
eval "$table_cmd"

sleep 30


#########################


# check that table_annovar completed

if [ ! -s "$annovar_multianno" ] ; then
	echo -e "\n $script_name ERROR: $annovar_multianno IS EMPTY \n" >&2
	exit 1
fi


#########################


# add identifier, remove some columns, and sort ANNOVAR table

# columns (original columns shifted by 1 by added identifier):
# 1 - mutation identifier
# 2-4 - Chr,Start,End
# 5-6 - Ref,Alt
# 7 - cytoBand
# 8-12 - refGene

# backslashes in awk to prevent variable expansion and retain quotes
bash_cmd="
cat $annovar_multianno \
| awk -F $'\t' 'BEGIN {OFS=FS} {print \$1 \":\" \$2 \"-\" \$3 , \$0}' \
| cut -f $annovar_keep_cols \
| sed 's/Chr:Start-End/#ID/g' \
| LC_ALL=C sort -k1,1 \
> $annovar_out_fixed
"
echo -e "\n CMD: $bash_cmd \n"
eval "$bash_cmd"

sleep 30


#########################


# check that table_annovar completed

if [ ! -s "$annovar_out_fixed" ] ; then
	echo -e "\n $script_name ERROR: $annovar_out_fixed IS EMPTY \n" >&2
	exit 1
fi


#########################


# convert regions table to a format for joining

# new header (with ID and sample name)
regions_table_cols=$(cat $regions_table | grep 'start.*end' | head -1)
echo -e "#ID\tSAMPLE\t${regions_table_cols}" > "$regions_table_fixed"

# table contents
cat $regions_table | grep -v 'start.*end' \
| awk -v sample="$sample" -F $'\t' 'BEGIN {OFS=FS} {print $1":"$2"-"$3,sample,$0}' \
| LC_ALL=C sort -k1,1 \
>> "$regions_table_fixed"

sleep 30


#########################


# check that table_annovar completed

if [ ! -s "$regions_table_fixed" ] ; then
	echo -e "\n $script_name ERROR: $regions_table_fixed IS EMPTY \n" >&2
	exit 1
fi


#########################


# merge variant info from VCF with annotations from ANNOVAR

join_cmd="
LC_ALL=C join -a1 -t $'\t' $regions_table_fixed $annovar_out_fixed > $annovar_combined
"
echo -e "\n CMD: $join_cmd \n"
eval "$join_cmd"

sleep 30


#########################


# check that join completed

if [ ! -s "$annovar_combined" ] ; then
	echo -e "\n $script_name ERROR: $annovar_combined IS EMPTY \n" >&2
	exit 1
fi


#########################


# clean up

rm -fv "$annovar_input"
rm -fv "$annovar_out_fixed"


#########################


# combine annotations for all samples

combine_all_cmd="
cat ${annovar_dir}/*.combined.txt | LC_ALL=C sort -k1,1 -k2,2 | uniq > ${annovar_dir}.txt
"
echo -e "\n CMD: $combine_all_cmd \n"
eval "$combine_all_cmd"


#########################



# end
