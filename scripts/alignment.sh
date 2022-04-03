#!/bin/sh
# The default partition is the 'general' partition
#SBATCH --partition=general
# The default Quality of Service is the 'short' QoS (maximum run time: 4 hours)
#SBATCH --qos=long
# The default run (wall-clock) time is 1 minute
#SBATCH --time=32:00:00
# The default number of parallel tasks per job is 1
#SBATCH --ntasks=1
# The default number of CPUs per task is 1 (note: CPUs are always allocated per 2)
# Request 1 CPU per active thread of your program (assume 1 unless you specifically set this)
#SBATCH --cpus-per-task=32
# The default memory per node is 1024 megabytes (1GB) (for multiple tasks, specify --mem-per-cpu instead)
#SBATCH --mem=80000
#SBATCH --mail-type=BEGIN,END
#SBATCH --chdir="/tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline"
#SBATCH --output=out_alignment_HG01258.out
#SBATCH --error=err_alignment_HG01258.out
#SBATCH --tmp=160G

export GENOME_ID="HG01258"
#source ~/.bashrc
# add channels to conda configuration
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# create the environment and install dependencies
conda create -y -n deepvariant_whatshap
conda activate deepvariant_whatshap
conda install -y whatshap==1.0 samtools==1.10 pbmm2 pbsv bcftools sniffles


#cd /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline
# index and align reference.
#samtools faidx reference/GRCh38_no_alt_analysis_set.fasta #this is done only once
#pbmm2 index reference/GRCh38_no_alt_analysis_set.fasta reference/ref.mmi --log-level INFO #this is done only once

#Need to align reads to refrence
#fofn way:

cd /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline/reads/${GENOME_ID}
touch myfiles.fofn
ls *.bam > myfiles.fofn
sed -i -e "s#^#reads/${GENOME_ID}/#" myfiles.fofn #add a prefix

cd /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline
samtools merge reads/${GENOME_ID}/${GENOME_ID}.bam -b reads/${GENOME_ID}/myfiles.fofn -f
pbmm2 align reference/ref.mmi reads/${GENOME_ID}/${GENOME_ID}.bam reads/${GENOME_ID}/aligned_reads.bam -m 2G -j 32 --preset CCS --sort --log-level INFO 