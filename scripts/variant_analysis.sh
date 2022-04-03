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
#SBATCH --output=out_analysis.out
#SBATCH --error=err_analysis.out
#SBATCH --tmp=50G


export GENOME_ID="NA12878"
#make consistent naming scheme
#both PBSV and Deepvariant have chrX in it, switch, do awk
#extract chromosomes from variant files

gzip -d deepvariant2/${GENOME_ID}/deepvariant_variants.vcf.gz
#REMOVE chr from chr14 etc...
awk '{gsub(/^chr/,""); print}' deepvariant2/${GENOME_ID}/deepvariant_variants.vcf > deepvariant2/${GENOME_ID}/deepvariant_variants_alt.vcf
awk '{gsub(/^chr/,""); print}' deepvariant2/${GENOME_ID}/deepvariant_variants_nosnp_filtered.vcf.recode.vcf > deepvariant2/${GENOME_ID}/deepvariant_variants_nosnp.vcf 
awk '{gsub(/^chr/,""); print}' pbsv/${GENOME_ID}/pbsv_variants.vcf > pbsv/${GENOME_ID}/pbsv_variants_alt.vcf
awk '{gsub(/^chr/,""); print}' sniffles/${GENOME_ID}/sniffles_variants.vcf > sniffles/${GENOME_ID}/sniffles_variants_alt.vcf

awk '{gsub(/^chr/,""); print}' variants/NA12878/NA12878.sorted.vcf > variants/NA12878/NA12878.sorted_alt.vcf

#Add chr to vcf
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' no_chr.vcf > with_chr.vcf
#https://www.biostars.org/p/98582/



#filter chr 6 and 14
vcftools --vcf deepvariant2/${GENOME_ID}/deepvariant_variants_alt.vcf --stdout --chr 6 --chr 14 --recode > deepvariant2/${GENOME_ID}/deepvariant_variants_filtered.vcf
vcftools --vcf pbsv/${GENOME_ID}/pbsv_variants_alt.vcf --stdout --chr 6 --chr 14 --recode > pbsv/${GENOME_ID}/pbsv_variants_alt_filtered.vcf 
vcftools --vcf sniffles/${GENOME_ID}/sniffles_variants_alt.vcf --stdout --chr 6 --chr 14 --recode > sniffles/${GENOME_ID}/sniffles_variants_filtered.vcf

#sort
vcf-sort pbsv/${GENOME_ID}/pbsv_variants_alt_filtered.vcf > pbsv/${GENOME_ID}/pbsv_variants_s.vcf
vcf-sort sniffles/${GENOME_ID}/sniffles_variants_filtered.vcf > sniffles/${GENOME_ID}/sniffles_variants_s.vcf

#zips
bgzip -c pbsv/${GENOME_ID}/pbsv_variants_s.vcf > pbsv/${GENOME_ID}/pbsv_variants.vcf.gz
bgzip -c sniffles/${GENOME_ID}/sniffles_variants_s.vcf > sniffles/${GENOME_ID}/sniffles_variants.vcf.gz
bgzip -c deepvariant2/${GENOME_ID}/deepvariant_variants_filtered.vcf > deepvariant2/${GENOME_ID}/deepvariant_variants.vcf.gz
bgzip -c deepvariant2/${GENOME_ID}/deepvariant_variants_nosnp.vcf > deepvariant2/${GENOME_ID}/deepvariant_variants_nosnp.vcf.gz

#tabix
tabix -p vcf pbsv/${GENOME_ID}/pbsv_variants.vcf.gz
tabix -p vcf sniffles/${GENOME_ID}/sniffles_variants.vcf.gz
tabix -p vcf deepvariant2/${GENOME_ID}/deepvariant_variants.vcf.gz
tabix -p vcf deepvariant2/${GENOME_ID}/deepvariant_variants_nosnp.vcf.gz


#vcf compare the 3 files
vcf-compare -w 2 pbsv/${GENOME_ID}/pbsv_variants.vcf.gz \
	deepvariant2/${GENOME_ID}/deepvariant_variants.vcf.gz \
	sniffles/${GENOME_ID}/sniffles_variants.vcf.gz > variants/${GENOME_ID}/compare_3_pipelines.out

#make a vcf merge file 

vcf-concat pbsv/${GENOME_ID}/pbsv_variants.vcf.gz \
	deepvariant2/${GENOME_ID}/deepvariant_variants.vcf.gz \
	sniffles/${GENOME_ID}/sniffles_variants.vcf.gz  | bgzip -c > variants/${GENOME_ID}/merged_variants.vcf.gz 

#compare vcf merge to some sort of gold standard IF WE HAVE IT
vcf-compare -w 2 variants/${GENOME_ID}/merged_variants.vcf.gz \ 
	variants/${GENOME_ID}/NA12878.sorted.vcf.gz > variants/${GENOME_ID}/compare_ground_truth.out
	

#process the reference
gzip -d variants/${GENOME_ID}/NA12878.sorted.vcf.gz -c > variants/${GENOME_ID}/NA12878.sorted.vcf
awk '{gsub(/^chr/,""); print}' variants/${GENOME_ID}/NA12878.sorted.vcf > variants/${GENOME_ID}/NA12878.sorted_alt.vcf
#vcftools --vcf variants/${GENOME_ID}/NA12878.sorted_alt.vcf --stdout --chr 6 --chr 14 --recode > variants/${GENOME_ID}/NA12878_filtered.vcf

bgzip -c variants/${GENOME_ID}/NA12878_filtered.vcf > variants/${GENOME_ID}/NA12878_filtered.vcf.gz
tabix -p vcf variants/${GENOME_ID}/NA12878_filtered.vcf.gz

vcf-compare -w 2 pbsv/${GENOME_ID}/pbsv_variants.vcf.gz \
	deepvariant2/${GENOME_ID}/deepvariant_variants_nosnp.vcf.gz \
	sniffles/${GENOME_ID}/sniffles_variants.vcf.gz \
	variants/${GENOME_ID}/NA12878_filtered.vcf.gz
	
	
gzip -d variants/${GENOME_ID}/chr14_1kg.vcf.gz
#Check, maybe we need to change chromosomes using awk
#head -n 50 variants/${GENOME_ID}/chr14_1kg.vcf
#rename chromosomes?

#remove snps
vcftools --vcf variants/${GENOME_ID}/chr14_1kg.vcf --stdout --keep-only-indels --chr 14 --recode > variants/${GENOME_ID}/chr14_1kg_filtered.vcf

#zips
bgzip -c variants/${GENOME_ID}/chr14_1kg.vcf > variants/${GENOME_ID}/chr14_1kg_filtered.vcf.gz

#tabix
tabix -p vcf variants/${GENOME_ID}/chr14_1kg_filtered.vcf.gz

vcf-compare -w 25 pbsv/${GENOME_ID}/pbsv_variants.vcf.gz deepvariant2/${GENOME_ID}/deepvariant_variants_nosnp.vcf.gz sniffles/${GENOME_ID}/sniffles_variants.vcf.gz variants/${GENOME_ID}/chr14_1kg_filtered.vcf.gz

#PBSV can separate DEL, INS, BND, INS.DUP
#Sniffles: <DUP>, <IND>, <DEL>

#separate <DEL> and <DUP>
cat <(grep '^#' pbsv_variants_alt.vcf) <(grep '<DEL>\|<DUP>' pbsv_variants_alt.vcf) > cnvs.vcf
cat <(grep '^#' pbsv_variants_alt.vcf) <(grep 'INS' pbsv_variants_alt.vcf) > insertions.vcf 

cat <(grep '^#' sniffles_variants_filtered.vcf) <(grep '<INS>' sniffles_variants_filtered.vcf) > insertions_s.vcf #for snffles
cat <(grep '^#' sniffles_variants_filtered.vcf) <(grep 'INV' sniffles_variants_filtered.vcf) > inversions.vcf #for snffles

#sniffles_variants_after_sort.vcf

export GENOME_ID="NA12878"
export VARIANT_TYPE="insertions"

vcf-sort pbsv/${GENOME_ID}/${VARIANT_TYPE}.vcf > pbsv/${GENOME_ID}/${VARIANT_TYPE}_s.vcf
vcf-sort sniffles/${GENOME_ID}/${VARIANT_TYPE}.vcf > sniffles/${GENOME_ID}/${VARIANT_TYPE}_s.vcf
#vcf-sort deepvariant2/${GENOME_ID}/${VARIANT_TYPE}.vcf > deepvariant2/${GENOME_ID}/${VARIANT_TYPE}_s.vcf

bgzip -c pbsv/${GENOME_ID}/${VARIANT_TYPE}_s.vcf > pbsv/${GENOME_ID}/${VARIANT_TYPE}.vcf.gz
bgzip -c sniffles/${GENOME_ID}/${VARIANT_TYPE}_s.vcf > sniffles/${GENOME_ID}/${VARIANT_TYPE}.vcf.gz
bgzip -c deepvariant2/${GENOME_ID}/${VARIANT_TYPE}_10bp.vcf > deepvariant2/${GENOME_ID}/${VARIANT_TYPE}_10bp.vcf.gz

tabix -p vcf pbsv/${GENOME_ID}/${VARIANT_TYPE}.vcf.gz
tabix -p vcf sniffles/${GENOME_ID}/${VARIANT_TYPE}.vcf.gz
tabix -p vcf deepvariant2/${GENOME_ID}/${VARIANT_TYPE}_10bp.vcf.gz

vcf-compare -w 25 pbsv/${GENOME_ID}/${VARIANT_TYPE}.vcf.gz sniffles/${GENOME_ID}/${VARIANT_TYPE}.vcf.gz deepvariant2/${GENOME_ID}/${VARIANT_TYPE}_10bp.vcf.gz

#bcftools view -v indels | head -n 200

#Deepvariant
#insertions	
bcftools view --types indels deepvariant2/NA12878/deepvariant_variants_nosnp.vcf.gz |
  bcftools norm -m - |
  bcftools filter --include 'strlen(REF)<strlen(ALT)' > deepvariant2/NA12878/insertions.vcf
  
#indels
bcftools view --types indels deepvariant2/NA12878/deepvariant_variants_nosnp_chr.vcf |
  bcftools norm -m - |
  bcftools filter --include 'strlen(REF) > 9 || strlen(ALT) > 9' > deepvariant2/NA12878/deepvariant_variants_nosnp_chr_10p.vcf

tabix -p vcf deepvariant2/${GENOME_ID}/deepvariant_variants_nosnp.vcf.gz


conda activate deepvariant_whatshap
#Regions
# 1 directory
# 2 filename

function create_filtered_vcfs {
	#bedtools intersect -a $1/${GENOME_ID}/$2.vcf -b benchmark/benchmark_repeats.bed -wa -header > $1/${GENOME_ID}/repeats.vcf
	#bedtools intersect -a $1/${GENOME_ID}/$2.vcf -b benchmark/benchmark_gap.bed -wa -header > $1/${GENOME_ID}/gaps.vcf
	#bedtools intersect -a $1/${GENOME_ID}/$2.vcf -b benchmark/benchmark_centromeres.bed -wa -header > $1/${GENOME_ID}/centromeres.vcf
	vcftools --vcf $1/${GENOME_ID}/$2.vcf --bed benchmark/benchmark_repeats.bed --recode --stdout > $1/${GENOME_ID}/repeats.vcf
	vcftools --vcf $1/${GENOME_ID}/$2.vcf --bed benchmark/benchmark_gap.bed --recode --stdout > $1/${GENOME_ID}/gaps.vcf
	vcftools --vcf $1/${GENOME_ID}/$2.vcf --bed benchmark/benchmark_centromeres.bed --recode --stdout > $1/${GENOME_ID}/centromeres.vcf
	
	bcftools norm -d none $1/${GENOME_ID}/repeats.vcf -o $1/${GENOME_ID}/repeats_normed.vcf
	bcftools norm -d none $1/${GENOME_ID}/gaps.vcf -o $1/${GENOME_ID}/gaps_normed.vcf
	bcftools norm -d none $1/${GENOME_ID}/centromeres.vcf -o $1/${GENOME_ID}/centromeres_normed.vcf

	#bedtools intersect -a $1/${GENOME_ID}/$2.vcf -b benchmark/gencode.v38.annotation.sorted.gff3 -wa -header > $1/${GENOME_ID}/coding.vcf
	#bedtools subtract -a $1/${GENOME_ID}/$2.vcf -b benchmark/gencode.v38.annotation.sorted.gff3 -wa -header > $1/${GENOME_ID}/noncoding.vcf

	vcftools --vcf $1/${GENOME_ID}/$2.vcf --bed benchmark/gencode.bed --recode --stdout > $1/${GENOME_ID}/coding.vcf
	vcftools --vcf $1/${GENOME_ID}/$2.vcf --exclude-bed benchmark/gencode.bed --recode --stdout > $1/${GENOME_ID}/noncoding.vcf
	bcftools norm -d none $1/${GENOME_ID}/coding.vcf -o $1/${GENOME_ID}/coding_normed.vcf
	bcftools norm -d none $1/${GENOME_ID}/noncoding.vcf -o $1/${GENOME_ID}/noncoding_normed.vcf
	
	#Create a stats TXT
	touch $1/${GENOME_ID}/region_stats 
	echo $2 > region_stats #filename
	
	echo 'repeats' >> $1/${GENOME_ID}/region_stats 
	grep -v "^#" $1/${GENOME_ID}/repeats_normed.vcf|wc -l >> $1/${GENOME_ID}/region_stats
	
	echo 'gaps' >> $1/${GENOME_ID}/region_stats
	grep -v "^#" $1/${GENOME_ID}/gaps_normed.vcf|wc -l >> $1/${GENOME_ID}/region_stats
	
	echo 'centromeres' >> $1/${GENOME_ID}/region_stats
	grep -v "^#" $1/${GENOME_ID}/centromeres_normed.vcf|wc -l >> $1/${GENOME_ID}/region_stats
	
	echo 'coding' >> $1/${GENOME_ID}/region_stats
	grep -v "^#" $1/${GENOME_ID}/coding_normed.vcf|wc -l >> $1/${GENOME_ID}/region_stats

	echo 'non-coding' >> $1/${GENOME_ID}/region_stats
	grep -v "^#" $1/${GENOME_ID}/noncoding_normed.vcf|wc -l >> $1/${GENOME_ID}/region_stats
}

create_filtered_vcfs deepvariant2 deepvariant_variants_nosnp_chr
create_filtered_vcfs pbsv pbsv_variants_chr 
create_filtered_vcfs sniffles sniffles_variants_chr 

# 1 directory
# 2 genomeID
# 3 filename
function create_length_vcfs {
	vcftools --vcf  $1/$2/$3.vcf --hist-indel-len --out $1/$2/lengths
}
#Transition into R for analysis

create_length_vcfs deepvariant2 ${GENOME_ID} deepvariant_variants_nosnp_chr_10p
create_length_vcfs pbsv ${GENOME_ID} pbsv_variants_chr 
create_length_vcfs sniffles ${GENOME_ID} sniffles_variants_chr 

vcftools --vcf  benchmark/NA12878_latest_benchmark_filtered.vcf --hist-indel-len --out benchmark/NA12878_latest_benchmark_lengths


function create_indels {
	#insertions
	bcftools view --types indels deepvariant2/${GENOME_ID}/$1.vcf.gz |
	bcftools norm -m - |
	bcftools filter --include 'strlen(REF)<strlen(ALT)' > deepvariant2/${GENOME_ID}/insertions.vcf

	cat <(grep '^#' pbsv/${GENOME_ID}/$2.vcf) <(grep 'INS' pbsv/${GENOME_ID}/$2.vcf) > pbsv/${GENOME_ID}/insertions.vcf
	
	cat <(grep '^#' sniffles/${GENOME_ID}/$3.vcf) <(grep '<INS>' sniffles/${GENOME_ID}/$3.vcf) > sniffles/${GENOME_ID}/insertions.vcf #for snffles
		
	#deletions
	bcftools view --types indels deepvariant2/${GENOME_ID}/$1.vcf.gz |
	bcftools norm -m - |
	bcftools filter --include 'strlen(REF)>strlen(ALT)' > deepvariant2/${GENOME_ID}/deletions.vcf
	
	cat <(grep '^#' pbsv/${GENOME_ID}/$2.vcf) <(grep 'DEL' pbsv/${GENOME_ID}/$2.vcf) > pbsv/${GENOME_ID}/deletions.vcf
	
	cat <(grep '^#' sniffles/${GENOME_ID}/$3.vcf) <(grep '<DEL>' sniffles/${GENOME_ID}/$3.vcf) > sniffles/${GENOME_ID}/deletions.vcf #for snffles
	
}


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
	
	touch $1/${GENOME_ID}/DP_stats 
	echo $2 > DP_stats #filename
	
	echo '10-' > $1/${GENOME_ID}/DP_stats 
	grep -v "^#" $1/${GENOME_ID}/$2.10dp.vcf|wc -l >> $1/${GENOME_ID}/DP_stats
	
	echo '10-19' >> $1/${GENOME_ID}/DP_stats
	grep -v "^#" $1/${GENOME_ID}/$2.20dp.vcf |wc -l >> $1/${GENOME_ID}/DP_stats
	
	echo '20-29' >> $1/${GENOME_ID}/DP_stats
	grep -v "^#" $1/${GENOME_ID}/$2.30dp.vcf |wc -l >> $1/${GENOME_ID}/DP_stats
	
	echo '30-49' >> $1/${GENOME_ID}/DP_stats
	grep -v "^#" $1/${GENOME_ID}/$2.50dp.vcf |wc -l >> $1/${GENOME_ID}/DP_stats

	echo '50-99' >> $1/${GENOME_ID}/DP_stats
	grep -v "^#" $1/${GENOME_ID}/$2.100dp.vcf |wc -l >> $1/${GENOME_ID}/DP_stats
	
	echo '100-199' >> $1/${GENOME_ID}/DP_stats
	grep -v "^#" $1/${GENOME_ID}/$2.200dp.vcf |wc -l >> $1/${GENOME_ID}/DP_stats
	
	echo '200+' >> $1/${GENOME_ID}/DP_stats
	grep -v "^#" $1/${GENOME_ID}/$2.large_dp.vcf |wc -l >> $1/${GENOME_ID}/DP_stats
}

create_dps deepvariant2 deepvariant_variants_nosnp_chr_10p


#1 is filename with path
#2 is end filename
function create_dps_comps {	
	touch cmp
	echo $1.10dp.vcf > cmp
	echo 'benchmark/NA12878_latest_benchmark_filtered.vcf' >> cmp
	
	./SURVIVOR/Debug/SURVIVOR merge cmp 1000 1 1 1 0 10 $1.10dp.overlap.vcf
	perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' $1.10dp.overlap.vcf | sed -e 's/\(.\)/\1 /g' > variants/NA12878/$2.10dp_overlap.txt
	#-----------------------------------------------------------------
	touch cmp
	echo $1.20dp.vcf > cmp
	echo 'benchmark/NA12878_latest_benchmark_filtered.vcf' >> cmp
	
	./SURVIVOR/Debug/SURVIVOR merge cmp 1000 1 1 1 0 10 $1.20dp.overlap.vcf
	perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' $1.20dp.overlap.vcf | sed -e 's/\(.\)/\1 /g' > variants/NA12878/$2.20dp_overlap.txt
	#-----------------------------------------------------------------
	touch cmp
	echo $1.30dp.vcf > cmp
	echo 'benchmark/NA12878_latest_benchmark_filtered.vcf' >> cmp
	
	./SURVIVOR/Debug/SURVIVOR merge cmp 1000 1 1 1 0 10 $1.30dp.overlap.vcf
	perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' $1.30dp.overlap.vcf | sed -e 's/\(.\)/\1 /g' > variants/NA12878/$2.30dp_overlap.txt
	#-----------------------------------------------------------------
	touch cmp
	echo $1.50dp.vcf > cmp
	echo 'benchmark/NA12878_latest_benchmark_filtered.vcf' >> cmp
	
	./SURVIVOR/Debug/SURVIVOR merge cmp 1000 1 1 1 0 10 $1.50dp.overlap.vcf
	perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' $1.50dp.overlap.vcf | sed -e 's/\(.\)/\1 /g' > variants/NA12878/$2.50dp_overlap.txt
	#-----------------------------------------------------------------
	touch cmp
	echo $1.100dp.vcf > cmp
	echo 'benchmark/NA12878_latest_benchmark_filtered.vcf' >> cmp

	./SURVIVOR/Debug/SURVIVOR merge cmp 1000 1 1 1 0 10 $1.100dp.overlap.vcf
	perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' $1.100dp.overlap.vcf | sed -e 's/\(.\)/\1 /g' > variants/NA12878/$2.100dp_overlap.txt
	#-----------------------------------------------------------------
	touch cmp
	echo $1.200dp.vcf > cmp
	echo 'benchmark/NA12878_latest_benchmark_filtered.vcf' >> cmp

	./SURVIVOR/Debug/SURVIVOR merge cmp 1000 1 1 1 0 10 $1.200dp.overlap.vcf
	perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' $1.200dp.overlap.vcf | sed -e 's/\(.\)/\1 /g' > variants/NA12878/$2.200dp_overlap.txt
	#-----------------------------------------------------------------

	touch cmp
	echo $1.large_dp.vcf > cmp
	echo 'benchmark/NA12878_latest_benchmark_filtered.vcf' >> cmp

	./SURVIVOR/Debug/SURVIVOR merge cmp 1000 1 1 1 0 10 $1.large_dp.overlap.vcf
	perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' $1.large_dp.overlap.vcf | sed -e 's/\(.\)/\1 /g' > variants/NA12878/$2.large_dp_overlap.txt
	#-----------------------------------------------------------------
}


#Basic stats, useful for getting DP
bcftools -stats -s- FILENAME.vcf 


#### 
####
#### GTF analysis
####
####

###Merging variant files using SURVIVOR
https://www.biostars.org/p/459269/

./SURVIVOR/Debug/SURVIVOR filter sniffles/NA12878/sniffles_variants_filtered_chr.vcf benchmark/benchmark_repeats.bed 20 -1 0.01 5 sniffles/NA12878/repeats.filtered.vcf

touch files_to_merge
echo 'pbsv/NA12878/pbsv_variants_chr.vcf' > files_to_merge
echo 'deepvariant2/NA12878/deepvariant_variants_nosnp_chr_10p.vcf' >> files_to_merge
echo 'sniffles/NA12878/sniffles_variants_filtered_chr_2.vcf' >> files_to_merge

touch files_to_cmpre
echo 'variants/NA12878/survivor_total_2.vcf' > files_to_cmpre
echo 'benchmark/HG001_benchmark_SV_indels.vcf' >> files_to_cmpre

touch two_benchmarks
echo 'variants/NA12878/NA12878_filtered_no_snp.vcf' > two_benchmarks
echo 'benchmark/HG001_benchmark_SV.vcf' >> two_benchmarks


#UNION of SVs 10 bp long
./SURVIVOR/Debug/SURVIVOR merge files_to_merge 1000 1 1 1 0 10 variants/NA12878/survivor_total_1.vcf
#INTERSECTION (2)
./SURVIVOR/Debug/SURVIVOR merge files_to_merge 1000 1 1 1 0 10 variants/NA12878/survivor_total_2.vcf

## compare to 'ground truth'
./SURVIVOR/Debug/SURVIVOR merge files_to_cmpre 1000 1 1 1 0 10 variants/NA12878/survivor_overlap_1_2.vcf

./SURVIVOR/Debug/SURVIVOR merge two_benchmarks 1000 1 1 1 0 10 variants/NA12878/survivor_benchmark_overlap.vcf

perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' variants/NA12878/survivor_benchmark_overlap.vcf | sed -e 's/\(.\)/\1 /g' > variants/NA12878/survivor_benchmark_overlap.txt

#IN R
t=read.table("FILENAME.txt",header=F)
library(VennDiagram)
venn.diagram(list(N15=which(t[,1]==1), T15=which(t[,2]==1)), fill = c("gray", "orange") , alpha = c(0.5, 0.5), cex = 2, lty =2, filename = "my_sample")

#bgzip -c variants/NA12878/NA12878.sorted.test.vcf > variants/NA12878/NA12878.sorted.test.vcf.gz

awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' variants/NA12878/NA12878_filtered.vcf > variants/NA12878/NA12878_filtered_chr.vcf

#vcftools --vcf variants/NA12878/NA12878.sorted.vcf  --stdout --chr 6 --chr 14 --recode > variants/NA12878/NA12878.chr6_14.vcf 

vcf-sort variants/NA12878/survivor_total_1.vcf > variants/NA12878/survivor_total_1_s.vcf

bgzip -c variants/NA12878/NA12878_filtered_chr.vcf > variants/NA12878/NA12878_filtered_chr.vcf.gz
bgzip -c variants/NA12878/survivor_total_1_s.vcf > variants/NA12878/survivor_total_1_s.vcf.gz

#tabix
tabix -p vcf variants/NA12878/NA12878_filtered_chr.vcf.gz
tabix -p vcf benchmark/NA12878_phased_variants.vcf.gz
tabix -p vcf variants/NA12878/survivor_total_1_s.vcf.gz

#vcf-compare -w 2 variants/NA12878/NA12878_filtered_chr.vcf.gz variants/NA12878/survivor_total_1_s.vcf.gz

#filter all variants shorter than 1000
bcftools view --types indels variants/NA12878/survivor_total_1_s.vcf.gz |
  bcftools norm -m - |
  bcftools filter --include 'strlen(REF) > 1000 || strlen(ALT) > 1000' > variants/NA12878/survivor_total_1_length_1k.vcf

bgzip -d variants/NA12878/survivor_total_1_length_1k.vcf.gz > variants/NA12878/survivor_total_1_length_1k.vcf

vcftools --gzvcf variants/NA12878/NA12878_filtered_chr.vcf.gz --stdout --keep-only-indels --chr 14 --chr 6 --recode > variants/NA12878/NA12878_filtered_chr_no_snp.vcf

#vcf-compare -w 200 benchmark/NA12878_phased_variants.vcf.gz variants/NA12878/survivor_total_1_s.vcf.gz #
vcf-compare -w 200 benchmark/HG001_benchmark_SV.vcf.gz variants/NA12878/survivor_total_1_s.vcf.gz #8.8% overlap


#bgzip -d benchmark/NA12878_phased_variants.vcf.gz > benchmark/NA12878_phased_variants.vcf.gz
vcftools --gzvcf benchmark/NA12878_phased_variants.vcf.gz --stdout --chr 14 --chr 6 --recode > benchmark/NA12878_phased_variants_indels.vcf
bgzip -c benchmark/NA12878_phased_variants_indels.vcf > benchmark/NA12878_phased_variants_indels.vcf.gz
tabix -p vcf benchmark/NA12878_phased_variants_indels.vcf.gz

vcftools --gzvcf benchmark/HG001_benchmark_SV.vcf.gz --stdout --keep-only-indels > benchmark/HG001_benchmark_SV.vcf
bgzip -c benchmark/HG001_benchmark_SV.vcf > benchmark/HG001_benchmark_SV_indels.vcf.gz
tabix -p vcf benchmark/HG001_benchmark_SV_indels.vcf.gz

vcf-compare -w 200 benchmark/HG001_benchmark_SV_indels.vcf.gz variants/NA12878/survivor_total_1_s.vcf.gz #8.4% overlap

vcf-compare -w 200 benchmark/NA12878_phased_variants_indels.vcf.gz variants/NA12878/survivor_total_1_s.vcf.gz #14.7% overlap

#with SNVs
vcf-compare -w 200 benchmark/NA12878_phased_variants.vcf.gz variants/NA12878/survivor_total_1_s.vcf.gz #47.4% overlap, with SNVs

