#!/bin/sh
# The default partition is the 'general' partition
#SBATCH --partition=general
# The default Quality of Service is the 'short' QoS (maximum run time: 4 hours)
#SBATCH --qos=long
# The default run (wall-clock) time is 1 minute
#SBATCH --time=1:00:00
# The default number of parallel tasks per job is 1
#SBATCH --ntasks=1
# The default number of CPUs per task is 1 (note: CPUs are always allocated per 2)
# Request 1 CPU per active thread of your program (assume 1 unless you specifically set this)
#SBATCH --cpus-per-task=8
# The default memory per node is 1024 megabytes (1GB) (for multiple tasks, specify --mem-per-cpu instead)
#SBATCH --mem=10000
#SBATCH --mail-type=BEGIN,END
#SBATCH --chdir="/tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline"
#SBATCH --output=out_survivor.out
#SBATCH --error=err_survivor.out
#SBATCH --tmp=10G

cd /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline

export GENOME_ID="HG003"

conda activate deepvariant_whatshap

mkdir variants/${GENOME_ID}

#UNION of SVs 10 bp long
./SURVIVOR/Debug/SURVIVOR merge files_to_merge 1000 1 1 1 0 10 variants/${GENOME_ID}/pipeline_union.vcf
#UNION of SVs 10 bp long
#./SURVIVOR/Debug/SURVIVOR merge files_to_merge 1000 1 1 1 0 10 variants/${GENOME_ID}/pipeline_union.vcf
#INTERSECTION (2)
./SURVIVOR/Debug/SURVIVOR merge files_to_merge 1000 2 1 1 0 10 variants/${GENOME_ID}/pipeline_intersection.vcf