for dir in "deepvariant2"/*
do
   cp $dir/lengths.indel.hist ./lengths/$dir.txt
done

for dir in "deepvariant2"/*
do
   cp $dir/region_stats_deepvariant2 ./region_stats/$dir.txt
done

for dir in "deepvariant2"/*
do
   cp $dir/DP_stats_deepvariant2 ./dp_stats/$dir.txt
done

for dir in "pbsv"/*
do
   cp $dir/lengths.indel.hist ./lengths/$dir.txt
done

for dir in "pbsv"/*
do
   cp $dir/region_stats_pbsv ./region_stats/$dir.txt
done

for dir in "pbsv"/*
do
   cp $dir/DP_stats_pbsv ./dp_stats/$dir.txt
done


for dir in "sniffles"/*
do
   cp $dir/lengths.indel.hist ./lengths/$dir.txt
done

for dir in "sniffles"/*
do
   cp $dir/region_stats_sniffles ./region_stats/$dir.txt
done


for dir in "sniffles"/*
do
	#run the grep file, make a new file
	cat $dir/sniffles_variants_filtered.vcf | grep -Po 'SVLEN=\K[^;]+' > $dir/all_lengths.txt
	#cat $dir/sniffles_variants_filtered.vcf | grep -Po 'SVLEN=\K[^;]+'
	#copy it into lengths_sniffles
	cp $dir/all_lengths.txt ./sniffles_lengths/$dir.txt
done

#here $dir is the absolute path, not the relative one, which is fucked up



#Filtering step for sniffles

bcftools query -e 'INFO/RE <5' -f '%LINE' deletions.vcf -> FILENAME
bcftools view -h deletions.vcf > headerfile
paste -d ""headerfile file2 > finalFile