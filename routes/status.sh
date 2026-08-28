#!/bin/bash


##
## Cluster/Slurm status (check that common modules can be loaded and executed)
##


# specify maximum runtime for sbatch job
# SBATCHTIME=15:00

# standard route header (validate args and print settings)
code_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
source "${code_dir}/scripts/route-header.sh" "$@"


#########################


# check that sbatch can be run
if ! sbatch --version >/dev/null 2>&1; then
	echo -e "\n $script_name ERROR: sbatch cannot be executed \n" >&2
	exit 1
fi

# check modules

echo "testing fastqc"
module purge
module add fastqc/0.11.7
if ! fastqc --version >/dev/null 2>&1; then
	echo -e "\n $script_name ERROR: fastqc cannot be executed at $(which fastqc) \n" >&2
	exit 1
fi

echo "testing fastq_screen"
module purge
module add fastq_screen/0.16.0
if ! fastq_screen --version >/dev/null 2>&1; then
	echo -e "\n $script_name ERROR: fastq_screen cannot be executed \n" >&2
	exit 1
fi
if [ ! -s "$(perldoc -l GD::Graph)" ] ; then
	echo -e "\n $script_name ERROR: perl GD::Graph module not found \n" >&2
	exit 1
fi

echo "testing trimmomatic"
module purge
module add java/1.8
module add trimmomatic/0.40
trimmomatic_jar=$(find "${TRIMMOMATIC_ROOT}" -maxdepth 1 -name "trimmomatic-*.jar" | head -1)
if [ ! -d "$TRIMMOMATIC_ROOT" ]; then
	echo -e "\n $script_name ERROR: trimmomatic root $TRIMMOMATIC_ROOT does not exist \n" >&2
	exit 1
fi
if [ ! -s "$trimmomatic_jar" ]; then
	echo -e "\n $script_name ERROR: trimmomatic jar $trimmomatic_jar does not exist \n" >&2
	exit 1
fi
if ! java -jar "$trimmomatic_jar" -version >/dev/null 2>&1; then
	echo -e "\n $script_name ERROR: trimmomatic cannot be executed at $trimmomatic_jar \n" >&2
	exit 1
fi

echo "testing star"
module purge
module add star/2.7.3a
# star/2.7.3a loads samtools/1.9
if ! STAR --version >/dev/null 2>&1; then
	echo -e "\n $script_name ERROR: STAR cannot be executed at $(which STAR) \n" >&2
	exit 1
fi
if ! samtools --version >/dev/null 2>&1; then
	echo -e "\n $script_name ERROR: samtools cannot be executed at $(which samtools) \n" >&2
	exit 1
fi

echo "testing bwa"
module purge
module add bwa/0.7.17
sambamba_bin="/gpfs/share/apps/sambamba/1.0.1/sambamba"
if ! bwa 2>&1 | grep -q "Version"; then
	echo -e "\n $script_name ERROR: bwa cannot be executed at $(which bwa) \n" >&2
	exit 1
fi
if ! $sambamba_bin --version >/dev/null 2>&1; then
	echo -e "\n $script_name ERROR: sambamba cannot be executed at $sambamba_bin \n" >&2
	exit 1
fi

echo "testing bowtie+samtools"
module purge
module add bowtie2/2.5.3
module add condaenvs/gpu/samtools1.23
sambamba_bin="/gpfs/share/apps/sambamba/1.0.1/sambamba"
if ! bowtie2 --version >/dev/null 2>&1; then
	echo -e "\n $script_name ERROR: bowtie2 cannot be executed at $(which bowtie2) \n" >&2
	exit 1
fi
if ! $sambamba_bin --version >/dev/null 2>&1; then
	echo -e "\n $script_name ERROR: sambamba cannot be executed at $sambamba_bin \n" >&2
	exit 1
fi
if ! samtools --version >/dev/null 2>&1; then
	echo -e "\n $script_name ERROR: samtools cannot be executed at $(which samtools) \n" >&2
	exit 1
fi

echo "testing samtools+bedtools"
module purge
module add condaenvs/gpu/samtools1.23
module add bedtools/2.30.0
if ! samtools --version >/dev/null 2>&1; then
	echo -e "\n $script_name ERROR: samtools cannot be executed at $(which samtools) \n" >&2
	exit 1
fi
if ! bcftools --version >/dev/null 2>&1; then
	echo -e "\n $script_name ERROR: bcftools cannot be executed at $(which bcftools) \n" >&2
	exit 1
fi
if ! bgzip --version >/dev/null 2>&1; then
	echo -e "\n $script_name ERROR: bgzip cannot be executed at $(which bgzip) \n" >&2
	exit 1
fi
if ! bedtools --version >/dev/null 2>&1; then
	echo -e "\n $script_name ERROR: bedtools cannot be executed at $(which bedtools) \n" >&2
	exit 1
fi


#########################


date



# end
