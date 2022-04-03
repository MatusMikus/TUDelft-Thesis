#Download file
curl -o LOCALNAME FTP_FILENAME

#Isolate 6, 14 chromosomes, remove snps, output compressed
#index 
tabix -p vcf FILENAME
tabix -p vcf NA12878_latest_benchmark.vcf.gz
tabix -p vcf sniffles/NA12878/sniffles_variants_min10bp.vcf.gz

bcftools view FILENAME --regions chr6,chr14 > FILENAME_AFTER
bcftools view NA12878_latest_benchmark.vcf.gz --regions chr6,chr14 -V snps -O z > NA12878_latest_benchmark_filtered.vcf.gz
bcftools view sniffles/NA12878/sniffles_variants_min10bp.vcf.gz --regions chr6,chr14 -V snps -O v > sniffles/NA12878/sniffles_variants_min10bp_chr.vcf

bgzip -d FILENAME.gz > FILENAME.vcf
bgzip -d NA12878_latest_benchmark_filtered.vcf.gz > NA12878_latest_benchmark_filtered.vcf

vcftools --vcf sniffles/NA12878/sniffles_variants_min10bp.vcf --stdout --chr 14 --chr 6 --recode > sniffles/NA12878/sniffles_variants_min10bp_recoded.vcf

#Remove chr from vcf
#awk '{gsub(/^chr/,""); print}' FILENAME.vcf > FILENAME.removed.vcf

#Add chr to vcf
#awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' no_chr.vcf > with_chr.vcf

#UNION of SVs 10 bp long
#./SURVIVOR/Debug/SURVIVOR merge files_to_merge 1000 1 1 1 0 10 variants/NA12878/survivor_total_1.vcf

#compare files total 2 is intersection
#total 1 is union
touch files_to_cmpre
echo 'variants/NA12878/survivor_total_2.vcf' > files_to_cmpre
echo 'benchmark/NA12878/NA12878_latest_benchmark_filtered.vcf' >> files_to_cmpre

## compare to 'ground truth'
./SURVIVOR/Debug/SURVIVOR merge files_to_cmpre 1000 1 1 1 0 10 variants/NA12878/survivor_latest_overlap.vcf

perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' variants/NA12878/survivor_latest_overlap.vcf | sed -e 's/\(.\)/\1 /g' > variants/NA12878/survivor_benchmark_overlap_latest.txt

touch sniffles_cmp
echo 'sniffles/NA12878/sniffles_variants_min10bp.vcf' > sniffles_cmp
echo 'benchmark/NA12878_latest_benchmark_filtered.vcf' >> sniffles_cmp

./SURVIVOR/Debug/SURVIVOR merge sniffles_cmp 1000 1 1 1 0 10 variants/NA12878/sniffles_10_latest_overlap.vcf

perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' variants/NA12878/sniffles_10_latest_overlap.vcf | sed -e 's/\(.\)/\1 /g' > variants/NA12878/sniffles_10_latest_overlap.txt


#STRINGENT SET CODE
touch stringent_set
echo 'pbsv/NA12878/pbsv_variants_chr.vcf' > stringent_set
echo 'deepvariant2/NA12878/deepvariant_variants_nosnp_chr_10p.vcf' >> stringent_set
echo 'sniffles/NA12878/sniffles_variants_filtered_chr_2.vcf' >> stringent_set

#INTERSECTION (2)

./SURVIVOR/Debug/SURVIVOR merge stringent_set 1000 2 1 1 0 10 variants/NA12878/survivor_2_callers_overlap.vcf
./SURVIVOR/Debug/SURVIVOR filter variants/NA12878/survivor_2_callers_overlap.vcf NA -1 -1 -1 5 variants/NA12878/survivor_2_callers_overlap_filtered.vcf

touch files_to_cmpre
echo 'variants/NA12878/survivor_2_callers_overlap.vcf' > files_to_cmpre
echo 'benchmark/NA12878_latest_benchmark_filtered.vcf' >> files_to_cmpre

./SURVIVOR/Debug/SURVIVOR merge files_to_cmpre 1000 1 1 1 0 10 variants/NA12878/benchmark_2_callers_overlap.vcf

perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' variants/NA12878/benchmark_2_callers_overlap.vcf | sed -e 's/\(.\)/\1 /g' > variants/NA12878/benchmark_2_callers_overlap.txt

#DOWNLOAD THE FILE LOCALLY

#IN R
t=read.table("FILENAME.txt",header=F)
library(VennDiagram)
venn.diagram(list(N15=which(t[,1]==1), T15=which(t[,2]==1)), fill = c("gray", "orange") , alpha = c(0.5, 0.5), cex = 2, lty =2, filename = "my_sample")

function create_dps {	
	bcftools view -i 'DP<10' $1.vcf > $1.10dp.vcf 
	touch cmp
	echo $1.10dp.vcf > cmp
	echo 'benchmark/NA12878_latest_benchmark_filtered.vcf' >> cmp
	
	./SURVIVOR/Debug/SURVIVOR merge cmp 1000 1 1 1 0 10 $1.10dp.overlap.vcf
	perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' $1.10dp.overlap.vcf | sed -e 's/\(.\)/\1 /g' > variants/NA12878/10dp_overlap.txt
	#-----------------------------------------------------------------
	bcftools view -i 'DP>=10 && DP<20' $1.vcf > $1.20dp.vcf 
	touch cmp
	echo $1.20dp.vcf > cmp
	echo 'benchmark/NA12878_latest_benchmark_filtered.vcf' >> cmp
	
	./SURVIVOR/Debug/SURVIVOR merge cmp 1000 1 1 1 0 10 $1.20dp.overlap.vcf
	perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' $1.20dp.overlap.vcf | sed -e 's/\(.\)/\1 /g' > variants/NA12878/20dp_overlap.txt
	#-----------------------------------------------------------------
	bcftools view -i 'DP>=20 && DP<30' $1.vcf > $1.30dp.vcf 
	touch cmp
	echo $1.30dp.vcf > cmp
	echo 'benchmark/NA12878_latest_benchmark_filtered.vcf' >> cmp
	
	./SURVIVOR/Debug/SURVIVOR merge cmp 1000 1 1 1 0 10 $1.30dp.overlap.vcf
	perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' $1.30dp.overlap.vcf | sed -e 's/\(.\)/\1 /g' > variants/NA12878/30dp_overlap.txt
	#-----------------------------------------------------------------
	bcftools view -i 'DP>=30 && DP<50' $1.vcf > $1.50dp.vcf 
	touch cmp
	echo $1.50dp.vcf > cmp
	echo 'benchmark/NA12878_latest_benchmark_filtered.vcf' >> cmp
	
	./SURVIVOR/Debug/SURVIVOR merge cmp 1000 1 1 1 0 10 $1.50dp.overlap.vcf
	perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' $1.50dp.overlap.vcf | sed -e 's/\(.\)/\1 /g' > variants/NA12878/50dp_overlap.txt
	#-----------------------------------------------------------------
	bcftools view -i 'DP>=50 && DP<100' $1.vcf > $1.100dp.vcf 
	touch cmp
	echo $1.100dp.vcf > cmp
	echo 'benchmark/NA12878_latest_benchmark_filtered.vcf' >> cmp

	./SURVIVOR/Debug/SURVIVOR merge cmp 1000 1 1 1 0 10 $1.100dp.overlap.vcf
	perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' $1.100dp.overlap.vcf | sed -e 's/\(.\)/\1 /g' > variants/NA12878/100dp_overlap.txt
	#-----------------------------------------------------------------
	bcftools view -i 'DP>=100 && DP<200' $1.vcf > $1.200dp.vcf
	touch cmp
	echo $1.200dp.vcf > cmp
	echo 'benchmark/NA12878_latest_benchmark_filtered.vcf' >> cmp

	./SURVIVOR/Debug/SURVIVOR merge cmp 1000 1 1 1 0 10 $1.200dp.overlap.vcf
	perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' $1.200dp.overlap.vcf | sed -e 's/\(.\)/\1 /g' > variants/NA12878/200dp_overlap.txt
	#-----------------------------------------------------------------

	bcftools view -i 'DP>=200' $1.vcf > $1.large_dp.vcf 
	touch cmp
	echo $1.large_dp.vcf > cmp
	echo 'benchmark/NA12878_latest_benchmark_filtered.vcf' >> cmp

	./SURVIVOR/Debug/SURVIVOR merge cmp 1000 1 1 1 0 10 $1.large_dp.overlap.vcf
	perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' $1.large_dp.overlap.vcf | sed -e 's/\(.\)/\1 /g' > variants/NA12878/large_dp_overlap.txt
	#-----------------------------------------------------------------
}




bcftools view --types indels benchmark/${GENOME_ID}/NA12878_latest_benchmark_filtered.vcf |
  bcftools norm -m - |
  bcftools filter --include 'strlen(REF) > 9 || strlen(ALT) > 9' > benchmark/${GENOME_ID}/NA12878_latest_benchmark_filtered_10bp.vcf