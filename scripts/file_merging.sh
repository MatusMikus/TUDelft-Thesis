#Run from full_genome_pipeline


#1: dir_name 
#2: variant_file_name
#this creates sniffles_all, sniffles_deletions_all etc...
function create_total {
	touch survivor_file
	touch survivor_file_del
	touch survivor_file_ins
	
	for dir in $1/* 
	do
		test -r $dir/$2.vcf && echo "$dir/$2.vcf" >> survivor_file
		test -r $dir/insertions.vcf && echo "$dir/insertions.vcf" >> survivor_file_ins
		test -r $dir/deletions.vcf && echo "$dir/deletions.vcf" >> survivor_file_del
	done
	
	cat survivor_file
	cat survivor_file_del
	cat survivor_file_ins
	./SURVIVOR/Debug/SURVIVOR merge survivor_file 1000 1 1 1 0 10 variants/"$1_all.vcf"
	echo "starting insertions"
	./SURVIVOR/Debug/SURVIVOR merge survivor_file_ins 1000 1 1 1 0 10 variants/"$1_insertions_all.vcf"
	echo "starting deletions" 
	./SURVIVOR/Debug/SURVIVOR merge survivor_file_del 1000 1 1 1 0 10 variants/"$1_deletions_all.vcf"
	rm survivor_file
	rm survivor_file_del
	rm survivor_file_ins
}

create_total pbsv pbsv_variants_filtered 
create_total deepvariant2 deepvariant_variants_filtered
create_total sniffles sniffles_variants_filtered

touch survivor_ins_all

echo 'variants/deepvariant2_insertions_all.vcf' > survivor_ins_all
echo 'variants/pbsv_insertions_all.vcf' >> survivor_ins_all
echo 'variants/sniffles_insertions_all.vcf' >> survivor_ins_all

./SURVIVOR/Debug/SURVIVOR merge survivor_ins_all 1000 1 0 0 0 10 variants/"insertions_all.vcf"


rm survivor_ins_all

touch survivor_del_all

echo 'variants/deepvariant2_deletions_all.vcf' > survivor_del_all
echo 'variants/pbsv_deletions_all.vcf' >> survivor_del_all
echo 'variants/sniffles_deletions_all.vcf' >> survivor_del_all

./SURVIVOR/Debug/SURVIVOR merge survivor_del_all 1000 1 0 0 0 10 variants/"deletions_all.vcf"

rm survivor_del_all

touch all_files

echo 'variants/pbsv_all.vcf' > all_files
echo 'variants/deepvariant2_all.vcf' >> all_files
echo 'variants/sniffles_all.vcf' >> all_files

./SURVIVOR/Debug/SURVIVOR merge all_files 1000 1 1 1 0 10 variants/all_variants.vcf

rm all_files

mkdir -p results	

awk '{if($3-$2 >= 10) print}' benchmark/GRCh38.nr_insertions_chr.bed > benchmark/GRCh38.nr_insertions_chr_10.bed
awk '{if($3-$2 >= 10) print}' benchmark/GRCh38.nr_deletions_chr.bed > benchmark/GRCh38.nr_deletions_chr_10.bed

function create_intersections { #insertions or deletions
	#these 4 will get the full hits in database
	bedtools intersect -u -f 0.9 -r -wa -a benchmark/GRCh38.nr_$1_chr_10.bed -b variants/deepvariant2_$1_all.vcf > results/dv_$1_90_u.bed
	bedtools intersect -u -f 0.9 -r -wa -a benchmark/GRCh38.nr_$1_chr_10.bed -b variants/sniffles_$1_all.vcf > results/sniffles_$1_90_u.bed
	bedtools intersect -u -f 0.9 -r -wa -a benchmark/GRCh38.nr_$1_chr_10.bed -b variants/pbsv_$1_all.vcf > results/pbsv_$1_90_u.bed
	bedtools intersect -u -f 0.9 -r -wa -a benchmark/GRCh38.nr_$1_chr_10.bed -b variants/$1_all.vcf > results/all_$1_90_u.bed
	
	#these 4 will get hits in results
	# bedtools intersect -f 0.9 -r -wb -a benchmark/GRCh38.nr_$1_chr.bed -b variants/deepvariant2_$1_all.vcf > results/dv_$1_90.bed
	# bedtools intersect -f 0.9 -r -wb -a benchmark/GRCh38.nr_$1_chr.bed -b variants/sniffles_$1_all.vcf > results/sniffles_$1_90.bed
	# bedtools intersect -f 0.9 -r -wb -a benchmark/GRCh38.nr_$1_chr.bed -b variants/pbsv_$1_all.vcf > results/pbsv_$1_90.bed
	# bedtools intersect -f 0.9 -r -wb -a benchmark/GRCh38.nr_$1_chr.bed -b variants/$1_all.vcf > results/all_$1_90.bed
	
	touch results/$1.txt
	echo 'deepvariant benchmark coverage' > results/$1.txt
	cat results/dv_$1_90_u.bed | wc -l | tr -d '\n' >> results/$1.txt
	echo -n '/' >> results/$1.txt
	cat benchmark/GRCh38.nr_$1_chr_10.bed | wc -l >> results/$1.txt
	# echo 'deepvariant percentage hit' >> results/$1.txt
	# cat results/dv_$1_90.bed | wc -l | tr -d '\n'  >> results/$1.txt
	# echo -n '/' >> results/$1.txt
	# grep -v "^#" variants/deepvariant2_$1_all.vcf | wc -l >> results/$1.txt 
	
	echo 'pbsv benchmark coverage' >> results/$1.txt
	cat results/pbsv_$1_90_u.bed | wc -l | tr -d '\n' >> results/$1.txt
	echo -n '/' >> results/$1.txt
	cat benchmark/GRCh38.nr_$1_chr_10.bed | wc -l >> results/$1.txt
	# echo 'pbsv percentage hit' >>	 results/$1.txt
	# cat results/pbsv_$1_90.bed | wc -l | tr -d '\n'  >> results/$1.txt
	# echo -n '/' >> results/$1.txt
	# grep -v "^#" variants/pbsv_$1_all.vcf | wc -l >> results/$1.txt
	
	echo 'sniffles benchmark coverage' >> results/$1.txt
	cat results/sniffles_$1_90_u.bed | wc -l | tr -d '\n' >> results/$1.txt
	echo -n '/' >> results/$1.txt
	cat benchmark/GRCh38.nr_$1_chr_10.bed | wc -l >> results/$1.txt
	# echo 'sniffles percentage hit' >> results/$1.txt
	# cat results/sniffles_$1_90.bed | wc -l | tr -d '\n'  >> results/$1.txt
	# echo -n '/' >> results/$1.txt
	# grep -v "^#" variants/sniffles_$1_all.vcf | wc -l >> results/$1.txt
	
	echo 'total benchmark coverage' >> results/$1.txt
	cat results/all_$1_90_u.bed | wc -l | tr -d '\n' >> results/$1.txt
	echo -n '/' >> results/$1.txt
	cat benchmark/GRCh38.nr_$1_chr_10.bed | wc -l >> results/$1.txt
	# echo 'total percentage hit' >> results/$1.txt
	# cat results/all_$1_90.bed | wc -l | tr -d '\n'  >> results/$1.txt
	# echo -n '/' >> results/$1.txt
	# grep -v "^#" variants/$1_all.vcf | wc -l >> results/$1.txt
}

#how many overlaps I have from the total - can have 5 overlaps on one region
create_intersections insertions
create_intersections deletions


#this is all in reference file
bedtools intersect -u -f 0.5 -a benchmark/GRCh38.nr_insertions_chr_10.bed -b variants/deepvariant2_insertions_all.vcf > results/test.bed
#write full B files - this shouldnt be more than results
bedtools intersect -f 0.99 -r -wb -a benchmark/GRCh38.nr_insertions_chr.bed -b variants/insertions_all.vcf > results/test_wb.bed

bedtools intersect -u -a benchmark/GRCh38.nr_deletions_chr.bed -b variants/deepvariant2_insertions_all.vcf > results/test_u.bed
bedtools intersect -a benchmark/GRCh38.nr_deletions_chr.bed -b variants/deepvariant2_all.vcf > results/test.bed

awk '{if($3-$2 >= 5) print}' benchmark/GRCh38.nr_insertions_chr.bed > benchmark/GRCh38.nr_insertions_chr_5.bed