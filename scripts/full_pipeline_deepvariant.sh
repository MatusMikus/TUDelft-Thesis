#!/bin/sh
# The default partition is the 'general' partition
#SBATCH --partition=general
# The default Quality of Service is the 'short' QoS (maximum run time: 4 hours)
#SBATCH --qos=long
# The default run (wall-clock) time is 1 minute
#SBATCH --time=16:00:00
# The default number of parallel tasks per job is 1
#SBATCH --ntasks=1
# The default number of CPUs per task is 1 (note: CPUs are always allocated per 2)
# Request 1 CPU per active thread of your program (assume 1 unless you specifically set this)
#SBATCH --cpus-per-task=8
# The default memory per node is 1024 megabytes (1GB) (for multiple tasks, specify --mem-per-cpu instead)
#SBATCH --mem=60000
#SBATCH --mail-type=BEGIN,END
#SBATCH --chdir="/tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline"
#SBATCH --output=out_deepvariant_HG01175.out
#SBATCH --error=err_deepvariant_HG01175.out
#SBATCH --tmp=220G

#curl -o HG002.bam https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_SequelII_CCS_11kb/HG002.SequelII.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam
ulimit -u 10000 # https://stackoverflow.com/questions/52026652/openblas-blas-thread-init-pthread-create-resource-temporarily-unavailable/54746150#54746150
export BIN_VERSION="1.1.0"
export GENOME_ID="HG00673"
#ln -s /tmp/${USER}/full_genome_pipeline /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline
tmpdir="/tmp/${USER}/full_genome_pipeline"
mkdir --parents "${tmpdir}"

# Cleanup temporary folder
function clean_up {
  rm --recursive --force "$tmpdir" && echo "Clean up of $tmpdir completed successfully."
  find /tmp/ -user mmikus -exec rm -rf {} \; #remove all other tmp files
  exit
}

# Setup clean_up to run on exit
trap 'clean_up' EXIT

#cp -r /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline/reference "${tmpdir}"  #copy the reference
ln -s '/tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline/reference' "${tmpdir}/"
ln -s '/tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline/whatshap' "${tmpdir}/"
ln -s '/tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline/whatshap/HG00673' "${tmpdir}/"

#You may want to uncomment this to know what to clean up when the clean up fails:
echo "Temporary folder: $(hostname --short):${tmpdir}"

export SINGULARITY_TMPDIR="${tmpdir}/tmp"
export SINGULARITY_CACHEDIR="${tmpdir}/cache"
export SINGULARITY_LOCALCACHEDIR="${tmpdir}/local_cache"
export SINGULARITYENV_HGREF="reference/GRCh38_no_alt_analysis_set.fasta"
export SINGULARITY_DISABLE_CACHE=True
export HGREF='reference/GRCh38_no_alt_analysis_set.fasta'

cd "${tmpdir}"
mkdir -p deepvariant1/${GENOME_ID}
mkdir -p whatshap/${GENOME_ID}
mkdir -p deepvariant2/${GENOME_ID}
mkdir -p pbsv/${GENOME_ID}
mkdir -p sniffles/${GENOME_ID}
mkdir -p reads/${GENOME_ID}
mkdir -p happy

#Make local cache folders to avoid singularity caching in home directory
mkdir -p tmp
mkdir -p cache
mkdir -p local_cache

#source ~/.bashrc
# add channels to conda configuration
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# create the environment and install dependencies
conda create -y -n deepvariant_whatshap
conda activate deepvariant_whatshap
conda install -y whatshap==1.0 samtools==1.10 pbmm2 pbsv bcftools sniffles

# cd /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline

#index and align reference.
#samtools faidx reference/GRCh38_no_alt_analysis_set.fasta #this is done only once
#pbmm2 index reference/GRCh38_no_alt_analysis_set.fasta reference/ref.mmi --log-level INFO #this is done only once

#Need to align reads to refrence
#Single file way
samtools merge reads/${GENOME_ID}/${GENOME_ID}.bam -b reads/${GENOME_ID}/myfiles.fofn
pbmm2 align reference/ref.mmi reads/${GENOME_ID}/${GENOME_ID}.bam reads/${GENOME_ID}/aligned_reads.bam -m 2G --preset CCS --sort --log-level INFO 

#Multiple files way, requires *.fofn file
#pbmm2 align reference/ref.mmi reads/${GENOME_ID}/myfiles.fofn reads/${GENOME_ID}/aligned_reads.bam --preset CCS --sort --log-level INFO

cd "${tmpdir}"

singularity exec --pwd ${tmpdir} \
	--bind /usr/lib/locale/ \
	--home /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline \
	--workdir /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline \
	docker://google/deepvariant:${BIN_VERSION} \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type PACBIO \
	--intermediate_results_dir /tmp/${USER}/ \
    --ref reference/GRCh38_no_alt_analysis_set.fasta \
    --reads /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline/reads/${GENOME_ID}/aligned_reads.bam \
    --output_vcf ${tmpdir}/deepvariant1/${GENOME_ID}/all_reads_strict.vcf.gz \
	--output_gvcf ${tmpdir}/deepvariant1/${GENOME_ID}/all_reads_strict.gvcf.gz \
	--make_examples_extra_args vsc_min_fraction_indels=0.1,vsc_min_count_snps=4,vsc_min_fraction_snps=0.15 \
    --num_shards $(nproc) \
	--regions "chr6 chr14"
	
	#takes up to 15x BAM file size in memory

whatshap phase \
        --output whatshap/${GENOME_ID}/all_reads_strict.phased.vcf.gz \
        --reference reference/GRCh38_no_alt_analysis_set.fasta \
        deepvariant1/${GENOME_ID}/all_reads_strict.vcf.gz \
        /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline/reads/${GENOME_ID}/aligned_reads.bam
		
tabix -p vcf ${tmpdir}/whatshap/${GENOME_ID}/all_reads_strict.phased.vcf.gz

#update halfway through, after DV 1
cp -ur "$tmpdir" /tudelft.net/staff-umbrella/hifidatafeb2021/
echo "copying halfway through"


whatshap haplotag \
        --output whatshap/${GENOME_ID}/all_reads_strict.haplotagged.bam \
        --reference reference/GRCh38_no_alt_analysis_set.fasta \
        whatshap/${GENOME_ID}/all_reads_strict.phased.vcf.gz \
        /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline/reads/${GENOME_ID}/aligned_reads.bam

samtools index whatshap/${GENOME_ID}/all_reads_strict.haplotagged.bam
echo "After samtools index"
	
singularity exec --pwd ${tmpdir} \
	--bind /usr/lib/locale/ \
	--home /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline \
	--workdir /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline \
	docker://google/deepvariant:${BIN_VERSION} \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type PACBIO \
	--intermediate_results_dir /tmp/${USER}/ \
    --ref reference/GRCh38_no_alt_analysis_set.fasta \
    --reads whatshap/${GENOME_ID}/all_reads_strict.haplotagged.bam \
	--use_hp_information \
    --output_vcf ${tmpdir}/deepvariant2/${GENOME_ID}/deepvariant_variants_strict.vcf.gz \
	--output_gvcf ${tmpdir}/deepvariant2/${GENOME_ID}/deepvariant_variants_strict.gvcf.gz \
	--make_examples_extra_args vsc_min_fraction_indels=0.1,vsc_min_count_snps=4,vsc_min_fraction_snps=0.15 \
    --num_shards $(nproc) \
	--regions "chr6 chr14"

#copy all created files back and delete the directory
cp -ur "$tmpdir" /tudelft.net/staff-umbrella/hifidatafeb2021/
rm -rf /tmp/${USER} 
find /tmp/ -user mmikus -exec rm -rf {} \; #remove all other tmp files