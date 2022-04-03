#!/bin/sh
# The default partition is the 'general' partition
#SBATCH --partition=general
# The default Quality of Service is the 'short' QoS (maximum run time: 4 hours)
#SBATCH --qos=long
# The default run (wall-clock) time is 1 minute
#SBATCH --time=4:00:00
# The default number of parallel tasks per job is 1
#SBATCH --ntasks=1
# The default number of CPUs per task is 1 (note: CPUs are always allocated per 2)
# Request 1 CPU per active thread of your program (assume 1 unless you specifically set this)
#SBATCH --cpus-per-task=8
# The default memory per node is 1024 megabytes (1GB) (for multiple tasks, specify --mem-per-cpu instead)
#SBATCH --mem=60000
#SBATCH --mail-type=BEGIN,END
#SBATCH --chdir="/tudelft.net/staff-umbrella/hifidatafeb2021/full_genome_pipeline"
#SBATCH --output=out_postprocessing.out
#SBATCH --error=err_postprocessing.out
#SBATCH --tmp=50G

#INPUT: 3 RAW result files with, consistent naming for genome
#OUTPUT: Length stats, region stats, DP stats, venn diagram for overlaps

export GENOME_ID="HG005"

# #conda doesn't support bgzip
# bgzip -c pbsv/${GENOME_ID}/pbsv_variants.vcf > pbsv/${GENOME_ID}/pbsv_variants.vcf.gz

conda create -y -n deepvariant_whatshap
conda activate deepvariant_whatshap
conda install -y whatshap==1.0 samtools==1.10 pbmm2 pbsv bcftools sniffles

# # remove sniffles breakpoints
# cat <(grep '^#' sniffles/${GENOME_ID}/sniffles_variants.vcf) <(grep -v 'N[[]\|[]]' sniffles/${GENOME_ID}/sniffles_variants.vcf) > sniffles/${GENOME_ID}/sniffles_variants_test.vcf
# vcf-sort sniffles/${GENOME_ID}/sniffles_variants_test.vcf > sniffles/${GENOME_ID}/sniffles_variants_test_s.vcf
# bgzip -c sniffles/${GENOME_ID}/sniffles_variants_test_s.vcf > sniffles/${GENOME_ID}/sniffles_variants_test.vcf.gz

# tabix -p vcf deepvariant2/${GENOME_ID}/deepvariant_variants_strict.vcf.gz
# tabix -p vcf pbsv/${GENOME_ID}/pbsv_variants.vcf.gz
# tabix -p vcf sniffles/${GENOME_ID}/sniffles_variants_test.vcf.gz

# bcftools view pbsv/${GENOME_ID}/pbsv_variants.vcf.gz --regions chr6,chr14 -V snps -O v > pbsv/${GENOME_ID}/pbsv_variants_filtered.vcf
# bcftools view deepvariant2/${GENOME_ID}/deepvariant_variants_strict.vcf.gz --regions chr6,chr14 -V snps -O v > deepvariant2/${GENOME_ID}/deepvariant_variants_filtered_chr.vcf
# bcftools view sniffles/${GENOME_ID}/sniffles_variants_test.vcf.gz --regions chr6,chr14 -V snps -O v > sniffles/${GENOME_ID}/sniffles_variants_filtered.vcf

# bcftools view --types indels deepvariant2/${GENOME_ID}/deepvariant_variants_filtered_chr.vcf |
  # bcftools norm -m - |
  # bcftools filter --include 'strlen(REF) > 9 || strlen(ALT) > 9' > deepvariant2/${GENOME_ID}/deepvariant_variants_filtered.vcf


#need to sort 

vcf-sort sniffles/${GENOME_ID}/sniffles_variants_test.vcf > sniffles/${GENOME_ID}/sniffles_variants_test_s.vcf

bgzip -c deepvariant2/${GENOME_ID}/deepvariant_variants_filtered.vcf > deepvariant2/${GENOME_ID}/deepvariant_variants_filtered.vcf.gz
bgzip -c pbsv/${GENOME_ID}/pbsv_variants_filtered.vcf > pbsv/${GENOME_ID}/pbsv_variants_filtered.vcf.gz
bgzip -c sniffles/${GENOME_ID}/sniffles_variants_filtered.vcf > sniffles/${GENOME_ID}/sniffles_variants_filtered.vcf.gz

tabix -p vcf deepvariant2/${GENOME_ID}/deepvariant_variants_filtered.vcf.gz
tabix -p vcf pbsv/${GENOME_ID}/pbsv_variants_filtered.vcf.gz
tabix -p vcf sniffles/${GENOME_ID}/sniffles_variants_filtered.vcf.gz



#filter out low read support
cp sniffles/${GENOME_ID}/sniffles_variants_filtered.vcf sniffles/${GENOME_ID}/sniffles_variants_filtered_backup.vcf
bcftools query -e 'INFO/RE <5' -f '%LINE' sniffles/${GENOME_ID}/sniffles_variants_filtered.vcf > sniffles/${GENOME_ID}/sniffles_variants_filtered_re_5.vcf
bcftools view -h sniffles/${GENOME_ID}/sniffles_variants_filtered.vcf > sniffles/${GENOME_ID}/headerfile
cat sniffles/${GENOME_ID}/headerfile sniffles/${GENOME_ID}/sniffles_variants_filtered_re_5.vcf > sniffles/${GENOME_ID}/sniffles_variants_filtered.vcf

bgzip -c sniffles/${GENOME_ID}/sniffles_variants_filtered.vcf > sniffles/${GENOME_ID}/sniffles_variants_filtered.vcf.gz


function create_filtered_vcfs {
	#bedtools intersect -a $1/${GENOME_ID}/$2.vcf -b benchmark/benchmark_repeats.bed -wa -header > $1/${GENOME_ID}/repeats.vcf
	#bedtools intersect -a $1/${GENOME_ID}/$2.vcf -b benchmark/benchmark_gap.bed -wa -header > $1/${GENOME_ID}/gaps.vcf
	#bedtools intersect -a $1/${GENOME_ID}/$2.vcf -b benchmark/benchmark_centromeres.bed -wa -header > $1/${GENOME_ID}/centromeres.vcf
	vcftools --vcf $1/${GENOME_ID}/$2.vcf --bed benchmark/old/benchmark_repeats.bed --recode --stdout > $1/${GENOME_ID}/repeats.vcf
	vcftools --vcf $1/${GENOME_ID}/$2.vcf --bed benchmark/old/benchmark_gap.bed --recode --stdout > $1/${GENOME_ID}/gaps.vcf
	vcftools --vcf $1/${GENOME_ID}/$2.vcf --bed benchmark/old/benchmark_centromeres.bed --recode --stdout > $1/${GENOME_ID}/centromeres.vcf
	
	bcftools norm -d none $1/${GENOME_ID}/repeats.vcf -o $1/${GENOME_ID}/repeats_normed.vcf
	bcftools norm -d none $1/${GENOME_ID}/gaps.vcf -o $1/${GENOME_ID}/gaps_normed.vcf
	bcftools norm -d none $1/${GENOME_ID}/centromeres.vcf -o $1/${GENOME_ID}/centromeres_normed.vcf

	#bedtools intersect -a $1/${GENOME_ID}/$2.vcf -b benchmark/gencode.v38.annotation.sorted.gff3 -wa -header > $1/${GENOME_ID}/coding.vcf
	#bedtools subtract -a $1/${GENOME_ID}/$2.vcf -b benchmark/gencode.v38.annotation.sorted.gff3 -wa -header > $1/${GENOME_ID}/noncoding.vcf

	vcftools --vcf $1/${GENOME_ID}/$2.vcf --bed benchmark/old/gencode.bed --recode --stdout > $1/${GENOME_ID}/coding.vcf
	vcftools --vcf $1/${GENOME_ID}/$2.vcf --exclude-bed benchmark/old/gencode.bed --recode --stdout > $1/${GENOME_ID}/noncoding.vcf
	bcftools norm -d none $1/${GENOME_ID}/coding.vcf -o $1/${GENOME_ID}/coding_normed.vcf
	bcftools norm -d none $1/${GENOME_ID}/noncoding.vcf -o $1/${GENOME_ID}/noncoding_normed.vcf
	
	#Create a stats TXT
	touch $1/${GENOME_ID}/region_stats_$1
	echo $2 > $1/${GENOME_ID}/region_stats_$1 #filename
	
	echo 'repeats' >> $1/${GENOME_ID}/region_stats_$1 
	grep -v "^#" $1/${GENOME_ID}/repeats_normed.vcf|wc -l >> $1/${GENOME_ID}/region_stats_$1
	
	echo 'gaps' >> $1/${GENOME_ID}/region_stats_$1
	grep -v "^#" $1/${GENOME_ID}/gaps_normed.vcf|wc -l >> $1/${GENOME_ID}/region_stats_$1
	
	echo 'centromeres' >> $1/${GENOME_ID}/region_stats_$1
	grep -v "^#" $1/${GENOME_ID}/centromeres_normed.vcf|wc -l >> $1/${GENOME_ID}/region_stats_$1
	
	echo 'coding' >> $1/${GENOME_ID}/region_stats_$1
	grep -v "^#" $1/${GENOME_ID}/coding_normed.vcf|wc -l >> $1/${GENOME_ID}/region_stats_$1

	echo 'non-coding' >> $1/${GENOME_ID}/region_stats_$1
	grep -v "^#" $1/${GENOME_ID}/noncoding_normed.vcf|wc -l >> $1/${GENOME_ID}/region_stats_$1
	
	rm $1/${GENOME_ID}/repeats_normed.vcf
	rm $1/${GENOME_ID}/repeats.vcf
	rm $1/${GENOME_ID}/gaps_normed.vcf
	rm $1/${GENOME_ID}/gaps.vcf
	rm $1/${GENOME_ID}/centromeres_normed.vcf
	rm $1/${GENOME_ID}/centromeres.vcf
	rm $1/${GENOME_ID}/coding_normed.vcf
	rm $1/${GENOME_ID}/coding.vcf
	rm $1/${GENOME_ID}/noncoding_normed.vcf
	rm $1/${GENOME_ID}/noncoding.vcf
}

create_filtered_vcfs deepvariant2 deepvariant_variants_filtered
create_filtered_vcfs pbsv pbsv_variants_filtered
create_filtered_vcfs sniffles sniffles_variants_filtered 
#create_filtered_vcfs benchmark NA12878_latest_benchmark_filtered


# 1 directory=
# 2 filename
function create_length_vcfs {
	vcftools --vcf  $1/${GENOME_ID}/$2.vcf --hist-indel-len --out $1/${GENOME_ID}/lengths
}

create_length_vcfs deepvariant2 deepvariant_variants_filtered
create_length_vcfs pbsv pbsv_variants_filtered
create_length_vcfs sniffles sniffles_variants_filtered 

function create_indels {
	#insertions
	#bcftools view --types indels deepvariant2/${GENOME_ID}/$1.vcf.gz |
	#bcftools norm -m - |
	#bcftools filter --include 'strlen(REF)<strlen(ALT)' > deepvariant2/${GENOME_ID}/insertions.vcf

	#cat <(grep '^#' pbsv/${GENOME_ID}/$2.vcf) <(grep 'INS' pbsv/${GENOME_ID}/$2.vcf) > pbsv/${GENOME_ID}/insertions.vcf
		
	cat <(grep '^#' sniffles/${GENOME_ID}/$3.vcf) <(grep 'SVTYPE=INS' sniffles/${GENOME_ID}/$3.vcf) > sniffles/${GENOME_ID}/insertions.vcf #for snffles
	#cat <(grep '^#' sniffles/${GENOME_ID}/$3.vcf) <(grep -v 'N[[]\|[]]\|<DEL>\|<DUP>\|<INV>' sniffles/${GENOME_ID}/$3.vcf) > sniffles/${GENOME_ID}/insertions.vcf #for snffles
		
	#deletions
	#bcftools view --types indels deepvariant2/${GENOME_ID}/$1.vcf.gz |
	#bcftools norm -m - |
	#bcftools filter --include 'strlen(REF)>strlen(ALT)' > deepvariant2/${GENOME_ID}/deletions.vcf
	
	#cat <(grep '^#' pbsv/${GENOME_ID}/$2.vcf) <(grep 'DEL' pbsv/${GENOME_ID}/$2.vcf) > pbsv/${GENOME_ID}/deletions.vcf

	cat <(grep '^#' sniffles/${GENOME_ID}/$3.vcf) <(grep '<DEL>\|SVTYPE=DEL' sniffles/${GENOME_ID}/$3.vcf) > sniffles/${GENOME_ID}/deletions.vcf #for snffles.
	#cat <(grep '^#' sniffles/${GENOME_ID}/$3.vcf) <(grep '<DEL>' sniffles/${GENOME_ID}/$3.vcf) > sniffles/${GENOME_ID}/deletions.vcf #for snffles. Also need to inlcude svtype=DEl
}

create_indels deepvariant_variants_filtered pbsv_variants_filtered sniffles_variants_filtered


# 1 directory
# 2 filename
function create_dps {	
	bcftools view -i 'DP<10' $1/${GENOME_ID}/$2.vcf > $1/${GENOME_ID}/$2.10dp.vcf 
	bcftools view -i 'DP>=10 && DP<20' $1/${GENOME_ID}/$2.vcf > $1/${GENOME_ID}/$2.20dp.vcf 
	bcftools view -i 'DP>=20 && DP<30' $1/${GENOME_ID}/$2.vcf > $1/${GENOME_ID}/$2.30dp.vcf 
	bcftools view -i 'DP>=30 && DP<50' $1/${GENOME_ID}/$2.vcf > $1/${GENOME_ID}/$2.50dp.vcf 
	bcftools view -i 'DP>=50 && DP<100' $1/${GENOME_ID}/$2.vcf > $1/${GENOME_ID}/$2.100dp.vcf 
	bcftools view -i 'DP>=100 && DP<200' $1/${GENOME_ID}/$2.vcf > $1/${GENOME_ID}/$2.200dp.vcf
	bcftools view -i 'DP>=200' $1/${GENOME_ID}/$2.vcf > $1/${GENOME_ID}/$2.large_dp.vcf 	
	
	touch $1/${GENOME_ID}/DP_stats_$1
	echo $2 > DP_stats_$1 #filename
	
	echo '10-' > $1/${GENOME_ID}/DP_stats_$1
	grep -v "^#" $1/${GENOME_ID}/$2.10dp.vcf|wc -l >> $1/${GENOME_ID}/DP_stats_$1
	
	echo '10-19' >> $1/${GENOME_ID}/DP_stats_$1
	grep -v "^#" $1/${GENOME_ID}/$2.20dp.vcf |wc -l >> $1/${GENOME_ID}/DP_stats_$1
	
	echo '20-29' >> $1/${GENOME_ID}/DP_stats_$1
	grep -v "^#" $1/${GENOME_ID}/$2.30dp.vcf |wc -l >> $1/${GENOME_ID}/DP_stats_$1
	
	echo '30-49' >> $1/${GENOME_ID}/DP_stats_$1
	grep -v "^#" $1/${GENOME_ID}/$2.50dp.vcf |wc -l >> $1/${GENOME_ID}/DP_stats_$1

	echo '50-99' >> $1/${GENOME_ID}/DP_stats_$1
	grep -v "^#" $1/${GENOME_ID}/$2.100dp.vcf |wc -l >> $1/${GENOME_ID}/DP_stats_$1
	
	echo '100-199' >> $1/${GENOME_ID}/DP_stats_$1
	grep -v "^#" $1/${GENOME_ID}/$2.200dp.vcf |wc -l >> $1/${GENOME_ID}/DP_stats_$1
	
	echo '200+' >> $1/${GENOME_ID}/DP_stats_$1
	grep -v "^#" $1/${GENOME_ID}/$2.large_dp.vcf |wc -l >> $1/${GENOME_ID}/DP_stats_$1
	
	rm $1/${GENOME_ID}/$2.10dp.vcf 
	rm $1/${GENOME_ID}/$2.20dp.vcf 
	rm $1/${GENOME_ID}/$2.30dp.vcf 
	rm $1/${GENOME_ID}/$2.50dp.vcf 
	rm $1/${GENOME_ID}/$2.100dp.vcf 
	rm $1/${GENOME_ID}/$2.200dp.vcf
	rm $1/${GENOME_ID}/$2.large_dp.vcf 	
}

create_dps deepvariant2 deepvariant_variants_filtered
create_dps pbsv pbsv_variants_filtered
#create_dps benchmark NA12878_latest_benchmark_filtered
#create_dps sniffles sniffles_variants_filtered

#merge into 1 total file

export GENOME_ID="HG00673"

mkdir variants/${GENOME_ID}

function create_survivor { 
	touch files_to_merge
	echo 'pbsv/'${GENOME_ID}'/pbsv_variants_filtered.vcf' > files_to_merge
	echo 'deepvariant2/'${GENOME_ID}'/deepvariant_variants_filtered.vcf' >> files_to_merge
	echo 'sniffles/'${GENOME_ID}'/sniffles_variants_filtered.vcf' >> files_to_merge
	
	#UNION of SVs 10 bp long
	./SURVIVOR/Debug/SURVIVOR merge files_to_merge 1000 1 1 1 0 10 variants/${GENOME_ID}/pipeline_union.vcf
	#UNION of SVs 10 bp long
	#./SURVIVOR/Debug/SURVIVOR merge files_to_merge 1000 1 1 1 0 10 variants/${GENOME_ID}/pipeline_union.vcf
	#INTERSECTION (2)
	./SURVIVOR/Debug/SURVIVOR merge files_to_merge 1000 2 1 1 0 10 variants/${GENOME_ID}/pipeline_intersection.vcf
}

create_survivor

touch files_to_merge
echo 'pbsv/'${GENOME_ID}'/pbsv_variants_filtered.vcf' > files_to_merge
echo 'deepvariant2/'${GENOME_ID}'/deepvariant_variants_filtered.vcf' >> files_to_merge
echo 'sniffles/'${GENOME_ID}'/sniffles_variants_filtered.vcf' >> files_to_merge

./SURVIVOR/Debug/SURVIVOR merge files_to_merge 1000 1 0 0 0 10 variants/${GENOME_ID}/pipeline_union_na.vcf
./SURVIVOR/Debug/SURVIVOR merge files_to_merge 1000 3 0 0 0 10 variants/${GENOME_ID}/pipeline_triple_na.vcf

bcftools isec -p test --nfiles=1 'pbsv/'${GENOME_ID}'/pbsv_variants_filtered.vcf.gz' 'deepvariant2/'${GENOME_ID}'/deepvariant_variants_filtered.vcf.gz' 'sniffles/'${GENOME_ID}'/sniffles_variants_filtered.vcf.gz'

intersection of all 3
take union of all 3, subtract all intersections - subtract a V b, b V c, c V a






# #UNION of SVs 10 bp long
# ./SURVIVOR/Debug/SURVIVOR merge files_to_merge 1000 1 1 1 0 10 variants/${GENOME_ID}/pipeline_union.vcf
# #UNION of SVs 10 bp long
# #./SURVIVOR/Debug/SURVIVOR merge files_to_merge 1000 1 1 1 0 10 variants/${GENOME_ID}/pipeline_union.vcf
# #INTERSECTION (2)
# ./SURVIVOR/Debug/SURVIVOR merge files_to_merge 1000 2 1 1 0 10 variants/${GENOME_ID}/pipeline_intersection.vcf

#./SURVIVOR/Debug/SURVIVOR merge files_to_merge 1000 3 1 1 0 10 variants/${GENOME_ID}/pipeline_intersection_all.vcf


#TODO: Cleanup
#See what files remain open, delete them
#reads
#whatshap
rm -r whatshap/${GENOME_ID}/all_reads_strict.haplotagged.bam
rm -r reads/${GENOME_ID}/aligned_reads.bam
rm -r reads/${GENOME_ID}/sniffles_reads_sorted.bam
rm -r reads/${GENOME_ID}/sniffles_reads.bam	

export GENOME_ID="NA12878"

cp sniffles/${GENOME_ID}/sniffles_variants_filtered.vcf sniffles/${GENOME_ID}/sniffles_variants_filtered_backup.vcf
bcftools query -e 'INFO/RE <5' -f '%LINE' sniffles/${GENOME_ID}/sniffles_variants_filtered.vcf > sniffles/${GENOME_ID}/sniffles_variants_filtered_re_5.vcf
bcftools view -h sniffles/${GENOME_ID}/sniffles_variants_filtered.vcf > sniffles/${GENOME_ID}/headerfile
cat sniffles/${GENOME_ID}/headerfile sniffles/${GENOME_ID}/sniffles_variants_filtered_re_5.vcf > sniffles/${GENOME_ID}/sniffles_variants_filtered.vcf

bgzip -c sniffles/${GENOME_ID}/sniffles_variants_filtered.vcf > sniffles/${GENOME_ID}/sniffles_variants_filtered.vcf.gz

create_filtered_vcfs sniffles sniffles_variants_filtered
create_length_vcfs sniffles sniffles_variants_filtered
create_indels deepvariant_variants_filtered pbsv_variants_filtered sniffles_variants_filtered
create_survivor