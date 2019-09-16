#!/bin/bash


##
## RNA-seq differential gene expression using DESeq2
##


# script filename
script_path="${BASH_SOURCE[0]}"
script_name=$(basename "$script_path")
route_name=${script_name/%.sh/}
echo -e "\n ========== ROUTE: $route_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 1 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir \n" >&2
	exit 1
fi

# standard comparison route arguments
proj_dir=$(readlink -f "$1")

# additional settings
code_dir=$(dirname $(dirname "$script_path"))

# display settings
echo
echo " * proj_dir: $proj_dir "
echo " * code_dir: $code_dir "
echo


#########################


# check that inputs exist

if [ ! -d "$proj_dir" ] ; then
	echo -e "\n $script_name ERROR: DIR $proj_dir DOES NOT EXIST \n" >&2
	exit 1
fi

gtf=$(bash "${code_dir}/scripts/get-set-setting.sh" "${proj_dir}/settings.txt" REF-GTF);

if [ ! -s "$gtf" ] || [ ! "$gtf" ] ; then
	echo -e "\n $script_name ERROR: GTF $gtf DOES NOT EXIST \n" >&2
	exit 1
fi

groups_table="${proj_dir}/samples.groups.csv"

if [ ! -s "$groups_table" ] ; then
	echo -e "\n $script_name ERROR: GROUP TABLE $groups_table DOES NOT EXIST \n" >&2
	exit 1
fi

num_samples=$(cat "$groups_table" | grep -v "#SAMPLE" | wc -l)

if [ "$num_samples" -lt 3 ] ; then
	echo -e "\n $script_name ERROR: $num_samples is too few samples \n" >&2
	exit 1
fi

num_groups=$(cat "$groups_table" | grep -v "#SAMPLE" | cut -d "," -f 2 | sort | uniq | wc -l)

gene_info_table="${proj_dir}/genes.featurecounts.txt"

if [ ! -s "$gene_info_table" ] ; then
	echo -e "\n $script_name ERROR: GENE INFO $gene_info_table DOES NOT EXIST \n" >&2
	exit 1
fi

strand=$(bash "${code_dir}/scripts/get-set-setting.sh" "${proj_dir}/settings.txt" EXP-STRAND);
counts_table="${proj_dir}/quant.featurecounts.counts.${strand}.txt"

if [ ! -s "$counts_table" ] ; then
	echo -e "\n $script_name ERROR: COUNTS TABLE $counts_table DOES NOT EXIST \n" >&2
	exit 1
fi

# unload all loaded modulefiles
module purge
module add default-environment


#########################


# settings and files

dge_dir="${proj_dir}/DGE-DESeq2-${strand}-${num_samples}samples-${num_groups}groups"

mkdir -p "$dge_dir"

input_groups_table="${dge_dir}/input.groups.csv"
input_counts_table="${dge_dir}/input.counts.txt"


#########################


# exit if output exits already

if [ -s "$input_groups_table" ] ; then
	echo -e "\n $script_name ERROR: TABLE $input_groups_table ALREADY EXISTS \n" >&2
	exit 1
fi

if [ -s "$input_counts_table" ] ; then
	echo -e "\n $script_name ERROR: TABLE $input_counts_table ALREADY EXISTS \n" >&2
	exit 1
fi


#########################


echo -e "\n ========== set up inputs ========== \n"

echo
echo " * counts table: $counts_table "
echo " * groups table: $groups_table "
echo

bash_cmd="rsync -t $counts_table $input_counts_table"
echo "CMD: $bash_cmd"
($bash_cmd)

bash_cmd="rsync -t $groups_table $input_groups_table"
echo "CMD: $bash_cmd"
($bash_cmd)

sleep 3


#########################


echo -e "\n ========== test R environment ========== \n"

# load relevant modules
module add r/3.6.1

echo
echo " * R: $(readlink -f $(which R)) "
echo " * R version: $(R --version | head -1) "
echo " * Rscript: $(readlink -f $(which Rscript)) "
echo " * Rscript version: $(Rscript --version 2>&1) "
echo

Rscript --vanilla "${code_dir}/scripts/test-package.R" optparse
Rscript --vanilla "${code_dir}/scripts/test-package.R" mnormt
Rscript --vanilla "${code_dir}/scripts/test-package.R" limma

sleep 5


#########################


echo -e "\n ========== start analysis ========== \n"

cd "$dge_dir" || exit 1

# launch the analysis R script
bash_cmd="Rscript --vanilla ${code_dir}/scripts/dge-deseq2.R $gtf $input_counts_table $input_groups_table"
echo "CMD: $bash_cmd"
($bash_cmd)


#########################


date



# end
