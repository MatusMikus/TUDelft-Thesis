#!/bin/sh
# The default partition is the 'general' partition
#SBATCH --partition=general
# The default Quality of Service is the 'short' QoS (maximum run time: 4 hours)
#SBATCH --qos=long
# The default run (wall-clock) time is 1 minute
#SBATCH --time=8:00:00
# The default number of parallel tasks per job is 1
#SBATCH --ntasks=1
# The default number of CPUs per task is 1 (note: CPUs are always allocated per 2)
# Request 1 CPU per active thread of your program (assume 1 unless you specifically set this)
#SBATCH --cpus-per-task=8
# The default memory per node is 1024 megabytes (1GB) (for multiple tasks, specify --mem-per-cpu instead)
#SBATCH --mem=700000
#SBATCH --mail-type=BEGIN,END
#SBATCH --chdir="/tudelft.net/staff-umbrella/hifidatafeb2021"
#SBATCH --output=out_pbsv_HG01123.out
#SBATCH --error=err_pbsv_HG01123.out
#SBATCH --tmp=130G

ulimit -u 10000 # https://stackoverflow.com/questions/52026652/openblas-blas-thread-init-pthread-create-resource-temporarily-unavailable/54746150#54746150
export BIN_VERSION="1.1.0"
export GENOME_ID="HG01123"

#ln -s /tmp/${USER}/full_genome_pipeline /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline

#tmpdir="/tmp/${USER}/full_genome_pipeline"
#mkdir --parents "${tmpdir}"

# Cleanup temporary folder
function clean_up {
  #rm --recursive --force "$tmpdir" && echo "Clean up of $tmpdir completed successfully."
  find /tmp/ -user mmikus -exec rm -rf {} \; #remove all other tmp files
  exit
}

# Setup clean_up to run on exit
trap 'clean_up' EXIT

#cp -r /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline/reference "${tmpdir}"  #copy the reference
#ln -s '/tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline/reference' "${tmpdir}/"

#You may want to uncomment this to know what to clean up when the clean up fails:
#echo "Temporary folder: $(hostname --short):${tmpdir}"

#export SINGULARITY_TMPDIR="${tmpdir}/tmp"
#export SINGULARITY_CACHEDIR="${tmpdir}/cache"
#export SINGULARITY_LOCALCACHEDIR="${tmpdir}/local_cache"
#export SINGULARITYENV_HGREF="reference/GRCh38_no_alt_analysis_set.fasta"
#export SINGULARITY_DISABLE_CACHE=True
#export HGREF='reference/GRCh38_no_alt_analysis_set.fasta'

#cd "${tmpdir}"

#Make local cache folders to avoid singularity caching in home directory
#mkdir -p tmp
#mkdir -p cache
#mkdir -p local_cache

cd /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline
mkdir -p pbsv/${GENOME_ID}

#source ~/.bashrc
# add channels to conda configuration
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# create the environment and install dependencies
conda create -y -n deepvariant_whatshap
conda activate deepvariant_whatshap
conda install -y whatshap==1.0 samtools==1.10 pbmm2 pbsv

# cd /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline/reads/HG00438/
# samtools merge test_reads.bam m64043_200710_174426.ccs.bam m64043_200711_235708.ccs.bam m64043_200714_124814.ccs.bam m64043_200713_062240.ccs.bam

#index reference
#samtools faidx reference/GRCh38_no_alt_analysis_set.fasta 

#align reads to reference
pbmm2 index reference/GRCh38_no_alt_analysis_set.fasta reference/ref.mmi --log-level INFO
pbmm2 align reference/ref.mmi reads/${GENOME_ID}/m64043_200711_235708.ccs.bam reads/HG00438/al_m64043_200711_235708.ccs.bam --sort --log-level INFO

#pbmm2 align reference/ref.mmi reads/HG00438/test_reads.bam reads/HG00438/aligned_reads.bam --sort --log-level INFO

#identify SV structures on chromosome 6 and chromosome 14
pbsv discover --tandem-repeats reference/human_GRCh38_no_alt_analysis_set.trf.bed \
	reads/${GENOME_ID}/aligned_reads.bam \
	reads/${GENOME_ID}/aligned_reads.svsig.gz

pbsv call --ccs -A 2 -O 2 -P 20 -m 10 reference/GRCh38_no_alt_analysis_set.fasta \
	/tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline/reads/${GENOME_ID}/aligned_reads.svsig.gz \
	/tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline/pbsv/${GENOME_ID}/pbsv_variants.vcf

find /tmp/ -user mmikus -exec rm -rf {} \; #remove all other tmp files


#cp -ur $"tmpdir" /tudelft.net/staff-umbrella/hifidatafeb2021/ #copy new/updated files back to my project folder
	
#cp -r /tmp/${USER}/full_genome_pipeline /tudelft.net/staff-umbrella/hifidatafeb2021/
#rmdir -f /tmp/${USER}

#Reads were aligned to GRCh38-NoALT (see Section 9.1 “Genome references”) using 
#"pbmm2 --sort --preset CCS -L 0.1 -c 0" for CCS and 
#"pbmm2 --sort --median-filter" for CLR. 

#Variants were called using PBSV v2.3.0 (SMRT Link v9.0, PacBio). 
#The PBSV workflow was executed separately per sample (not jointly) and per chromosome. First,
#SV signatures were discovered with "pbsv discover --tandem-repeats <SR.bed> -r <CHROM>"

#where SR.bed is the GRCh38 UCSC simpleRepeats track with regions within 200 bp merged ("bedtools merge -d 200"). 
#Next, variants were called with 
#"pbsv call --ccs -O 2 -P 20 -m 10" for CCS 
#O Ignore calls supported by < N reads in every sample. [2]
#P Ignore calls supported by < P% of reads in every sample. [20]
#m min SV Lentgh
#"pbsv call -m 10" for CLR. 

#For each sample, the per-chromosome calls were concatenated, sorted, and compressed with BCFtools 1.9 (117)