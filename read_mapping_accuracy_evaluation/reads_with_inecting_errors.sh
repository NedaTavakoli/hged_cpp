# Read without injecting any errors

# select random variant position form the first 6000 lines
cat data/variant_positions_snps_indels_chr22.txt | head -6000 >  variant_positions_snps_indels_chr22_head_6000.txt 
shuf -n 5000 variant_positions_snps_indels_chr22_head_6000.txt > random_pos.txt 


num_reads=5000
variant_positions=($(cut -d ',' -f2 random_pos.txt))

samples=($($bcftools query -l data/chr22_snps_indels.vcf.gz)) # array of samples, index from 0
count=0
for v in "${variant_positions[@]}"
do
    # echo "Location"$v
    sample=${samples[$count]} # random sample
    count=$(($count+1))
    # echo "Haplotype: "$sample
    s1=$($samtools faidx data/hs37d5.fa 22:${v}-$((${v}+100-1)) | $bcftools consensus -s ${sample} -H 2 data/chr22_snps_indels.vcf.gz |  sed '1d' | tr -d "[:space:]")
    # echo $s1
    a=($v $sample $s1) 
    printf '%s\n' "${a[*]}" >> pos.hap.reads.chr22_${num_reads}.txt  
done
cat pos.hap.reads.chr22.txt  | awk '{print $3}' > list.reads.chr22_${num_reads}.txt 


./vg map -t 36 -x chr22_snps_indels.xg -g chr22_snps_indels.gcsa -T list.reads.chr22_100.txt  > chr22_snps_indels.gam
./vg surject -x chr22_snps_indels.xg -b chr22_snps_indels.gam > chr22_snps_indels.bam
$samtools flagstat chr22_snps_indels.bam


./vg map -t 36 -x reduced_chr22_snps_indels_100_1.xg -g reduced_chr22_snps_indels_100_1.gcsa -T list.reads.chr22_100.txt > reduced_chr22_snps_indels_100_1.gam
./vg surject -x reduced_chr22_snps_indels_100_1.xg -b reduced_chr22_snps_indels_100_1.gam > reduced_chr22_snps_indels_100_1.bam
$samtools flagstat reduced_chr22_snps_indels_100_1.bam


GGGAAGTGGAGC

./vg map -t 36 -x reduced_chr22_snps_indels_100_1.xg -g reduced_chr22_snps_indels_100_1.gcsa -s TTGCATTGAGCCGAGATTGCGCCACTGCAGTCAAAACAGTCCGGCCTGGGCGACAGAGCGAGACTCCGTCTC > single_on_hap.gam
./vg surject -x reduced_chr22_snps_indels_100_1.xg -b single_on_hap.gam > single_on_hap.bam
$samtools flagstat  single_on_hap.bam

./vg map -t 36 -x chr22_snps_indels.xg -g chr22_snps_indels.gcsa -s TTGCATTGAGCCGAGATTGCGCCACTGCAGTCAAAACAGTCCGGCCTGGGCGACAGAGCGAGACTCCGTCTC > single_on_hap_comp.gam
./vg surject -x chr22_snps_indels.xg -b single_on_hap_comp.gam > single_on_hap_comp.bam
$samtools view -h chr22_snps_indels.bam > chr22_snps_indels.sam
cat chr22_snps_indels.sam
$samtools flagstat chr22_snps_indels.bam

