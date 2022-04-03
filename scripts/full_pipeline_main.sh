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
#SBATCH --mem=60000
#SBATCH --mail-type=BEGIN,END
#SBATCH --chdir="/tudelft.net/staff-umbrella/hifidatafeb2021"
#SBATCH --output=out_deepvariant.out
#SBATCH --error=err_deepvariant.out
#SBATCH --tmp=50G

ulimit -u 10000 # https://stackoverflow.com/questions/52026652/openblas-blas-thread-init-pthread-create-resource-temporarily-unavailable/54746150#54746150
export BIN_VERSION="1.1.0"

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

#You may want to uncomment this to know what to clean up when the clean up fails:
echo "Temporary folder: $(hostname --short):${tmpdir}"

export SINGULARITY_TMPDIR="${tmpdir}/tmp"
export SINGULARITY_CACHEDIR="${tmpdir}/cache"
export SINGULARITY_LOCALCACHEDIR="${tmpdir}/local_cache"
export SINGULARITYENV_HGREF="reference/GRCh38_no_alt_analysis_set.fasta"
export SINGULARITY_DISABLE_CACHE=True
export HGREF='reference/GRCh38_no_alt_analysis_set.fasta'

cd "${tmpdir}"
mkdir -p deepvariant1/HG00438
mkdir -p whatshap/HG00438
mkdir -p deepvariant2/HG00438
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
conda install -y whatshap==1.0 samtools==1.10 pbmm2 pbsv


#cd /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline
# index and align reference.
#samtools faidx reference/GRCh38_no_alt_analysis_set.fasta #this is done only once
#pbmm2 index reference/GRCh38_no_alt_analysis_set.fasta reference/ref.mmi --log-level INFO #this is done only once
#pbmm2 align reference/ref.mmi reads/HG00438/test_reads.bam reads/HG00438/aligned_reads.bam --sort --log-level INFO

cd "${tmpdir}"

# #Align reads
# # cd reads
# # for dir in reads
	# # do 
		# # for file in "$/dir/*.bam"
		# # do
		  # # pbmm2 align ../reference/ref.mmi "$file" "al_$file"
		# # done
	# # done
# # cd ..

# singularity exec --pwd ${tmpdir} \
	# --bind /usr/lib/locale/ \
	# --home /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline \
	# --workdir /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline \
	# docker://google/deepvariant:${BIN_VERSION} \
    # /opt/deepvariant/bin/run_deepvariant \
    # --model_type PACBIO \
	# --intermediate_results_dir /tmp/${USER}/ \
    # --ref reference/GRCh38_no_alt_analysis_set.fasta \
    # --reads /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline/reads/HG00438/aligned_reads.bam \
    # --output_vcf ${tmpdir}/deepvariant1/HG00438/all_reads.vcf.gz \
	# --output_gvcf ${tmpdir}/deepvariant1/HG00438/all_reads.gvcf.gz \
    # --num_shards $(nproc) \
	# --regions "chr6 chr14"
	
	# #takes up to 15x BAM file size in memory

whatshap phase \
        --output whatshap/HG00438/all_reads.phased.vcf.gz \
        --reference reference/GRCh38_no_alt_analysis_set.fasta \
        deepvariant1/HG00438/all_reads.vcf.gz \
        /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline/reads/HG00438/aligned_reads.bam
		
tabix -p vcf ${tmpdir}/whatshap/HG00438/all_reads.phased.vcf.gz

whatshap haplotag \
        --output whatshap/HG00438/all_reads.haplotagged.bam \
        --reference reference/GRCh38_no_alt_analysis_set.fasta \
        whatshap/HG00438/all_reads.phased.vcf.gz \
        /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline/reads/HG00438/aligned_reads.bam

samtools index whatshap/HG00438/all_reads.haplotagged.bam
	
singularity exec --pwd ${tmpdir} \
	--bind /usr/lib/locale/ \
	--home /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline \
	--workdir /tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline \
	docker://google/deepvariant:${BIN_VERSION} \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type PACBIO \
	--intermediate_results_dir /tmp/${USER}/ \
    --ref reference/GRCh38_no_alt_analysis_set.fasta \
    --reads whatshap/HG00438/all_reads.haplotagged.bam \
	--use_hp_information \
    --output_vcf ${tmpdir}/deepvariant2/HG00438/all_reads.vcf.gz \
	--output_gvcf ${tmpdir}/deepvariant2/HG00438/all_reads.gvcf.gz \
    --num_shards $(nproc) \
	--regions "chr6 chr14"

#copy all created files back and delete the directory
cp -ur "$tmpdir" /tudelft.net/staff-umbrella/hifidatafeb2021/
rm -rf /tmp/${USER} 
find /tmp/ -user mmikus -exec rm -rf {} \; #remove all other tmp files

# For multiple files
# declare -a trio=(HG002 HG003 HG004)
# for SAMPLE in "${trio[@]}"
# do
  # BAM=${SAMPLE}.bam

  # OUTPUT_VCF=${SAMPLE}.vcf.gz
  # OUTPUT_GVCF=${SAMPLE}.g.vcf.gz

  # time sudo docker run \
    # -v "${DIR}":"/data" \
    # google/deepvariant:${VERSION} \
    # /opt/deepvariant/bin/run_deepvariant \
    # --model_type=WES \
    # --ref="/data/hs37d5.fa" \
    # --reads="/data/${BAM}" \
    # --regions="/data/${CAPTURE_BED}" \
    # --output_vcf="/data/${OUTPUT_VCF}" \
    # --output_gvcf="/data/${OUTPUT_GVCF}" \
    # --num_shards=${N_SHARDS}
# done

#cp -r /tmp/${USER}/full_genome_pipeline /tudelft.net/staff-umbrella/hifidatafeb2021/
