#!/bin/sh
# The default partition is the 'general' partition
#SBATCH --partition=general
# The default Quality of Service is the 'short' QoS (maximum run time: 4 hours)
#SBATCH --qos=long
# The default run (wall-clock) time is 1 minute
#SBATCH --time=24:00:00
# The default number of parallel tasks per job is 1
#SBATCH --ntasks=1
# The default number of CPUs per task is 1 (note: CPUs are always allocated per 2)
# Request 1 CPU per active thread of your program (assume 1 unless you specifically set this)
#SBATCH --cpus-per-task=8
# The default memory per node is 1024 megabytes (1GB) (for multiple tasks, specify --mem-per-cpu instead)
#SBATCH --mem=80000
#SBATCH --mail-type=BEGIN,END
#SBATCH --chdir="/tudelft.net/staff-umbrella/hifidatafeb2021"
#SBATCH --output=out_sniffles_HG01123.out
#SBATCH --error=err_sniffles_HG01123.out
#SBATCH --tmp=150G

ulimit -u 10000 # https://stackoverflow.com/questions/52026652/openblas-blas-thread-init-pthread-create-resource-temporarily-unavailable/54746150#54746150
export BIN_VERSION="1.1.0"
export GENOME_ID="HG01123"

source ~/.bashrc  

#ln -s /tmp/${USER}/full_genome_pipeline /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline
# tmpdir="/tmp/${USER}/full_genome_pipeline"
# mkdir --parents "${tmpdir}"


# Cleanup temporary folder
# function clean_up {
  # rm --recursive --force "$tmpdir" && echo "Clean up of $tmpdir completed successfully."
  # exit
# }

# # Setup clean_up to run on exit
# trap 'clean_up' EXIT

#cd /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline/reads/HG00438/
#samtools calmd -b aligned_reads.bam ../../reference/GRCh38_no_alt_analysis_set.fasta > al_sniffles.bam
mkdir -p /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline/sniffles/${GENOME_ID}

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# create the environment and install dependencies
conda create -y -n deepvariant_whatshap
conda clean -ay
conda activate deepvariant_whatshap
conda install -y whatshap==1.0 samtools==1.10 pbmm2 pbsv bcftools sniffles

mkdir -p /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline/sniffles/${GENOME_ID}
cd /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline

#for sniffles	
samtools calmd -b reads/${GENOME_ID}/aligned_reads.bam reference/GRCh38_no_alt_analysis_set.fasta > reads/${GENOME_ID}/sniffles_reads.bam
samtools sort -@ 4 -m 10G reads/${GENOME_ID}/sniffles_reads.bam -o reads/${GENOME_ID}/sniffles_reads_sorted.bam

#custom sniffles
#cd /tudelft.net/staff-umbrella/hifidatafeb2021/Sniffles-master/bin/sniffles-core-1.0.12/
#./sniffles --min_support 1 -l 20 --genotype --skip_parameter_estimation -m ../../../full_genome_pipeline/reads/${GENOME_ID}/sniffles_reads.bam -v ../../../full_genome_pipeline/sniffles/${GENOME_ID}/sniffles_variants_after_sort.vcf

#./sniffles --min_support 1 --skip_parameter_estimation -m ../../../full_genome_pipeline/reads/HG00438/al_sniffles.bam -v ../../../full_genome_pipeline/sniffles/HG00438/variants.vcf

sniffles --min_support 2 -l 10 --genotype --skip_parameter_estimation -m reads/${GENOME_ID}/sniffles_reads_sorted.bam -v sniffles/${GENOME_ID}/sniffles_variants.vcf