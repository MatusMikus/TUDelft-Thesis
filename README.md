Used tools:
aws cli v2\
Bioconda: whatshap samtools pbmm2 pbsv bcftools sniffles vcftools\


Use `full_pipeline_prep_simple.sh` to download the reference and individual genomes.\
`alignment.sh` was used to align to the reference from multiple files.\
Use `full_pipeline_deepvariant.sh`, `full_pipeline_pbsv.sh`, `full_pipeline_sniffles.sh` for each individual algorithm, setting the GENOME_ID variable for individual genomes.\
the remaining scripts were used for creation of intermediate files and analysis.\