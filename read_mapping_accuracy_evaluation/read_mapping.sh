
# 16050075 16654125
ssh ntavakoli6@login-phoenix-3.pace.gatech.edu
qsub -I -q inferno -A GT-saluru8-CODA20 -l nodes=1:ppn=1,mem=300gb,walltime=10:00:00

# hive 
samtools=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/samtools-1.12/samtools
bcftools=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/bcftools-1.9/bcftools
tabix=/storage/home/hhive1/ntavakoli6/data/hged/software/htslib-1.12/tabix
bgzip=/storage/home/hhive1/ntavakoli6/data/hged/software/htslib-1.12/bgzip


samtools=/storage/home/hcoda1/6/ntavakoli6/p-saluru8-0/hged/software/samtools-1.12/samtools
bcftools=/storage/home/hcoda1/6/ntavakoli6/p-saluru8-0/hged/software/bcftools-1.9/bcftools
#***********************************************************
# get dependencies
$bcftools view -v 'snps,indels' -Oz chr22.vcf.gz > chr22_snps_indels.vcf.gz
$bcftools index chr22_snps_indels.vcf.gz # or $tabix -p vcf chr22_snps_indels.vcf.gz
$bcftools view -H -v 'snps,indels' chr22.vcf.gz | awk '{print $2}' > variant_positions_snps_indels_chr22.txt
$bcftools view -H -v 'snps,indels' chr22.vcf.gz | awk '{print $2"\t"$4"\t"$5}' > chr22_snps_indel_POS_REF_ALT.txt
$bcftools query -l chr22_snps_indels.vcf.gz > samples.txt

#***********************************************************
# get variant retained vcf file

python src/solution_analyzer.py graph_chr22_10k_100.txt \
    ILP_sol_vectors/ILP_sol_100_1.txt retained_variants_ed_ILP_100_1.txt    

# do only once
cat chr22_snps_indels.vcf | grep '^#' > header_chr22.vcf
cat chr22_snps_indels.vcf | grep  -vE '^#' > non_header_chr22.txt

# do for each config
cp header_chr22.vcf retained_variants_ed_ILP_100_1_chr22.vcf
awk 'NR == FNR {a[$0]; next } $2 in a {print $0} ' retained_variants_ed_ILP_100_1.txt  non_header_chr22.txt >> retained_variants_ed_ILP_100_1_chr22.vcf

# Generate vcf.gz file and its index file vcf.gz.tbi
cp retained_variants_ed_ILP_100_1_chr22.vcf retained_variants_ed_ILP_100_1_chr22_copy.vcf
$bgzip -c retained_variants_ed_ILP_100_1_chr22.vcf > retained_variants_ed_ILP_100_1_chr22.vcf.gz
$tabix -p vcf retained_variants_ed_ILP_100_1_chr22.vcf.gz
#***********************************************************

# READ MAPPING startrs from here

# construct list of random;y chosen haplotypes at range and location
# location: select a random number from a variant position list 

# Read without injecting any errors
num_reads=50
cat data/variant_positions_snps_indels_chr22.txt | head -${num_reads} > data/variant_positions_snps_indels_chr22_head_${num_reads}.txt
variant_positions=($(cut -d ',' -f2 data/variant_positions_snps_indels_chr22_head_${num_reads}.txt))
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
    printf '%s\n' "${a[*]}" >> pos.hap.reads.chr22.txt  
done
cat pos.hap.reads.chr22.txt  | awk '{print $3}' > list.reads.chr22_${num_reads}.txt 


#----------
# Map on complete graph

# construct graph from snps and indels
./vg construct -r data/hs37d5.fa -v chr22_snps_indels.vcf.gz -a -f -m 32 > graph.chr22.vg
#index the graph
./vg index -x graph.chr22.xg -g graph.chr22.gcsa graph.chr22.vg

-bash-4.2$ cat pos_substrings_chr22_10k_100.txt  | head -1
16050075 GGTGGGCCTAAGTGCCTCCTCTCGGGACTGGTATGGGGACGGTCATGCAA

# to map a single read
./vg map -s GGTGGGCCTAAGTGCCTCCTCTCGGGACTGGTATGGGGACGGTCATGCAA -x graph.chr22.xg -g graph.chr22.gcsa > read.chr22.gam
./vg surject -x graph.chr22.xg -b read.chr22.gam > read.chr22.gam.bam
$samtools flagstat read.chr22.gam.bam


# to map lists of reads
./vg map -T list.reads.chr22.txt -x graph.chr22.xg -g graph.chr22.gcsa > list.read.chr22.gam
./vg surject -x graph.chr22.xg -b list.read.chr22.gam > list.read.chr22.gam.bam
$samtools flagstat list.read.chr2.bam

#----------
# on reduced graph
./vg construct -r data/hs37d5.fa -v retained_variants_ed_ILP_100_1_chr22.vcf.gz -a -f -m 32 > reduced_100_1.graph.chr22.vg
# warning:[vg::Constructor] Unsupported IUPAC ambiguity codes found in 3; coercing to N.
# warning:[vg::Constructor] Lowercase characters found in hs37d5; coercing to uppercase.
# warning:[vg::Constructor] Unsupported IUPAC ambiguity codes found in hs37d5; coercing to N.

#index the graph
./vg index -x  reduced_100_1.graph.chr22.xg -g  reduced_100_1.graph.chr22.gcsa  reduced_100_1.graph.chr22.vg

# to map a single read
./vg map -s GGTGGGCCTAAGTGCCTCCTCTCGGGACTGGTATGGGGACGGTCATGCAA -x reduced_100_1.graph.chr22.xg -g reduced_100_1.graph.chr22.gcsa > reduced_100_1.read.chr22.gam
./vg surject -x reduced_100_1.graph.chr22.xg -b reduced_100_1.read.chr22.gam > reduced_100_1.read.chr22.gam.bam
$samtools flagstat reduced_100_1.read.chr22.gam.bam


# to map lists of reads
./vg map -T list.reads.chr22.txt -x reduced_100_1.graph.chr22.xg -g reduced_100_1.graph.chr22.gcsa > reduced_100_1.list.read.chr22.gam
./vg surject -x reduced_100_1.graph.chr22.xg -b reduced_100_1.list.read.chr22.gam > reduced_100_1.list.read.chr22.gam.bam
$samtools flagstat reduced_100_1.list.read.chr22.bam



# Using pruned version

# construct graph
# vg construct -r hg38.fa -v chr2.vcf.gz -R chr2 -C > chr2.vg
# The chr2.vg is 1.6G.

# xg index
# vg index -t 32 -x chr2.xg chr2.vg
# The chr2.xg is 3.4G.

# prune
# vg prune -r -p -t 32 chr2.vg > chr2.pruned.vg
# The chr2.pruned.vg is 308M.

# gcsa index
# vg index -t 32 -g chr2.gcsa chr2.pruned.vg
# The chr2.gcsa is 2.0G and chr2.gcsa.lcp is 1.2G.

# map
# vg map -t 36 -x chr2.xg -g chr2.gcsa -f 1.fq -f 2.fq > chr2.gam

##------

./vg construct -r data/hs37d5.fa -v data/chr22_snps_indels.vcf.gz -R 22 -C > chr22_snps_indels.vg
./vg index -t 32 -x chr22_snps_indels.xg chr22_snps_indels.vg
./vg prune -r -p -t 32 chr22_snps_indels.vg > chr22_snps_indels.pruned.vg
#Original graph chr22_snps_indels.vg: 4328928 nodes, 5514038 edges
# Built a temporary XG index
# Removed all paths
# Pruned complex regions: 4328928 nodes, 5163819 edges
# Removed small subgraphs: 3974884 nodes, 4892472 edges
# Restored graph: 4212289 nodes
# Serialized the graph: 4212289 nodes, 5167863 edges
./vg index -t 32 -g chr22_snps_indels.gcsa chr22_snps_indels.pruned.vg
./vg map -t 36 -x chr22_snps_indels.xg -g chr22_snps_indels.gcsa -T list.reads.chr22_50.txt  > chr22_snps_indels.gam
./vg surject -x chr22_snps_indels.xg -b chr22_snps_indels.gam > chr22_snps_indels.bam
$samtools flagstat chr22_snps_indels.bam



./vg construct -r data/hs37d5.fa -v retained_variants_ed_ILP_100_1_chr22.vcf.gz -R 22 -C > reduced_chr22_snps_indels_100_1.vg
./vg index -t 32 -x reduced_chr22_snps_indels_100_1.xg reduced_chr22_snps_indels_100_1.vg
./vg prune -r -p -t 32 reduced_chr22_snps_indels_100_1.vg > reduced_chr22_snps_indels_100_1.pruned.vg
# Original graph reduced_chr22_snps_indels_100_1.vg: 1604635 nodes, 1605249 edges
# Built a temporary XG index
# Removed all paths
# Pruned complex regions: 1604635 nodes, 1605211 edges
# Removed small subgraphs: 1604598 nodes, 1605185 edges
# Restored graph: 1604623 nodes
# Serialized the graph: 1604623 nodes, 1605214 edges
./vg index -t 32 -g reduced_chr22_snps_indels_100_1.gcsa reduced_chr22_snps_indels_100_1.pruned.vg
./vg map -t 36 -x reduced_chr22_snps_indels_100_1.xg -g reduced_chr22_snps_indels_100_1.gcsa -T list.reads.chr22_50.txt > reduced_chr22_snps_indels_100_1.gam
./vg surject -x reduced_chr22_snps_indels_100_1.xg -b reduced_chr22_snps_indels_100_1.gam > reduced_chr22_snps_indels_100_1.bam
$samtools flagstat reduced_chr22_snps_indels_100_1.bam




# just for test: single sequence

./vg map -t 36 -x reduced_chr22_snps_indels_100_1.xg -g reduced_chr22_snps_indels_100_1.gcsa -s AAATTTTCTTTCTTCTCGGGGGA > one.gam
./vg surject -x reduced_chr22_snps_indels_100_1.xg -b one.gam > one.bam
$samtools flagstat one.bam

