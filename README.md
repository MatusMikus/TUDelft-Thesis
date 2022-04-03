Used tools:
aws cli v2\
Bioconda: whatshap samtools pbmm2 pbsv bcftools sniffles vcftools\


Use `full_pipeline_prep_simple.sh` to download the reference and individual genomes.\
Use `alignment.sh` in conjunction with the above script to align reads that span multiple files.\
Use `full_pipeline_deepvariant.sh`, `full_pipeline_pbsv.sh`, `full_pipeline_sniffles.sh` for each individual algorithm, setting the GENOME_ID variable for individual genomes.\
the remaining scripts were used for creation of intermediate files and analysis.\
\

Data repositories:\
`https://github.com/genome-in-a-bottle/giab_data_indexes`\
`https://github.com/human-pangenomics/HPP_Year1_Data_Freeze_v1.0/`