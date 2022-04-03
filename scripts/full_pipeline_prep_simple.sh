#!/bin/sh

# You can control the resources and scheduling with '#SBATCH' settings
# (see 'man sbatch' for more information on setting these parameters)

# The default partition is the 'general' partition
#SBATCH --partition=general
# The default Quality of Service is the 'short' QoS (maximum run time: 4 hours)
#SBATCH --qos=long
# The default run (wall-clock) time is 1 minute
#SBATCH --time=4:00:00
# The default number of parallel tasks per job is 1
#SBATCH --ntasks=1
# Request 1 CPU per active thread of your program (assume 1 unless you specifically set this)
# The default number of CPUs per task is 1 (note: CPUs are always allocated per 2)
#SBATCH --cpus-per-task=1
# The default memory per node is 1024 megabytes (1GB) (for multiple tasks, specify --mem-per-cpu instead)
#SBATCH --mem=4096
# Set mail type to 'END' to receive a mail when the job finishes
# Do not enable mails when submitting large numbers (>20) of jobs at once
#SBATCH --mail-type=END
#SBATCH --chdir="/tudelft.net/staff-umbrella/hifidatafeb2021"
#SBATCH --output=full_pipeline_out.out
#SBATCH --error=full_pipeline_err.out
#SBATCH --tmp=10G

export SSSDIR=s3://human-pangenomics/working/HPRC
export GENOME_ID="HG00673"
mkdir /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline/reads/${GENOME_ID}
aws s3 cp ${SSSDIR}/${GENOME_ID}/raw_data/PacBio_HiFi/ /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline/reads/${GENOME_ID} --recursive

#https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=working/HPRC/HG01258/raw_data/PacBio_HiFi_Swap_Fixed/

# export NCBIDIR="ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/sv_datasets/nonredundant/"
# export DELS="deletions/GRCh38.nr_deletions.bed.gz"
# export DUPS="duplications/GRCh38.nr_duplications.bed.gz"
# export INS="insertions/GRCh38.nr_insertions.bed.gz"

# curl ${NCBIDIR}${DELS} -O  

# gunzip GRCh38.nr_deletions.bed.gz; 
# echo "track name=\"dbVar NR deletions\" description=\"non-redundant deletions from dbVar\"" > GRCh38.nr_deletions_ucsc.bed; 
# grep -v ^chrMT GRCh38.nr_deletions.bed >> GRCh38.nr_deletions_ucsc.bed

# gunzip GRCh38.nr_duplications.bed.gz

# gunzip GRCh38.nr_insertions.bed.gz

# grep '^chr6*|^chr14*' GRCh38.nr_insertions.bed > GRCh38.nr_insertions_chr.bed

mkdir -p full_genome_pipeline
cd full_genome_pipeline

#source ~/.bashrc  

# index reference
mkdir -p reference
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | gunzip > reference/GRCh38_no_alt_analysis_set.fasta
#  samtools faidx reference/GRCh38_no_alt_analysis_set.fasta

mkdir -p deepvariant1
mkdir -p whatshap
mkdir -p deepvariant2
#mkdir -p happy #unnecessary for now

#Make local cache folders to avoid singularity caching in home directory
mkdir -p tmp
mkdir -p cache
mkdir -p local_cache