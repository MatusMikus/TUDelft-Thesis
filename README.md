### Used tools\    
aws cli v2, singularity containers\
Algorithms: PBSV, Sniffles, DeepVariant\
Bioconda: whatshap samtools pbmm2 pbsv bcftools sniffles vcftools bcftools\
R for statistics/analysis

### Use\
Use `full_pipeline_prep_simple.sh` to download the reference and individual genomes.\
Use `alignment.sh` in conjunction with the above script to align reads that span multiple files.\
Use `full_pipeline_deepvariant.sh`, `full_pipeline_pbsv.sh`, `full_pipeline_sniffles.sh` for each individual algorithm, setting the GENOME_ID variable for individual genomes.\
the remaining scripts were used for creation of intermediate files and analysis.\
`analysis.R` was used to create statistics

### Data repositories\
`https://github.com/genome-in-a-bottle/giab_data_indexes`\
`https://github.com/human-pangenomics/HPP_Year1_Data_Freeze_v1.0/`