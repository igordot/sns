#!/bin/bash


##
## route header (sourced by route scripts)
##


# script filename (0 for actual script, 1 for sourced from)
script_path="${BASH_SOURCE[1]}"

# show route info
script_name=$(basename "$script_path")
route_name=${script_name/%.sh/}
echo -e "\n ========== ROUTE: $route_name ========== \n" >&2

# check if $code_dir is missing or doesn't point to a real directory
if [[ -z "${code_dir:-}" || ! -d "$code_dir" ]]; then
	echo -e "\n $script_name ERROR: code directory not found \n" >&2
	exit 1
fi

# check for the correct number of arguments
if [ ! $# == 2 ] ; then
	echo -e "\n $script_name ERROR: script requires two arguments \n" >&2
	echo -e "\n USAGE: $script_name project_dir sample_name \n" >&2
	exit 1
fi

# set standard route arguments
proj_dir=$(readlink -f "$1")
sample=$2

# reserve a thread for overhead
threads=$SLURM_CPUS_PER_TASK
threads=$(( threads - 1 ))

# print settings
echo
echo " * proj_dir: $proj_dir "
echo " * sample: $sample "
echo " * code_dir: $code_dir "
echo " * slurm threads: $SLURM_CPUS_PER_TASK "
echo " * command threads: $threads "
echo " * slurm nodename: $SLURMD_NODENAME "
echo " * hostname: $(hostname) "
echo " * time: $(date "+%Y-%m-%d %H:%M") "
echo

# default-environment module is needed to run sbatch (not loaded by default on all nodes)
module purge
module add default-environment
