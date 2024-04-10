#!/bin/bash

#M2305834 NGS pipeline 

# Setting up directories
mkdir -p ~/final_ngs/dnaseq
cd ~/final_ngs/dnaseq
mkdir data meta results logs
cd data
mkdir untrimmed_fastq trimmed_fastq
cd untrimmed_fastq

# Downloading data
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

# Change file format to gz
mv NGS0001.R1.fastq.qz NGS0001.R1.fastq.gz
mv NGS0001.R2.fastq.qz NGS0001.R2.fastq.gz

# Install tools
cd ~/
wget https://repo.anaconda.com/archive/Anaconda3-2022.10-Linux-x86_64.sh
chmod +x ./Anaconda3-2022.10-Linux-x86_64.sh
bash ./Anaconda3-2022.10-Linux-x86_64.sh -b -p $HOME/anaconda3
source $HOME/anaconda3/bin/activate
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install samtools
conda install bwa
conda install freebayes
conda install picard
conda install bedtools
conda install trimmomatic
conda install fastqc
conda install vcflib
# Add any other tools you need here

# View files / Uncompressing
cd ~/final_ngs/dnaseq/data/untrimmed_fastq
zcat NGS0001.R1.fastq.gz > NGS0001.R1.fastq
grep -B1 -A2 NNNNNN NGS0001.R1.fastq > bad_reads.txt
zcat NGS0001.R2.fastq.gz > NGS0001.R2.fastq
grep -B1 -A2 NNNNNN NGS0001.R2.fastq >> bad_reads.txt
mkdir ../other/
mv bad_reads.txt ../other/
grep NNNNNN NGS0001.R1.fastq | wc
grep NNNNNN NGS0001.R1.fastq | wc -l

# Loop script for bad reads trimming

ls -l
cd ~/final_ngs/dnaseq/data/untrimmed_fastq
ls *fastq.gz
filenames=`ls *.fastq.gz`
echo $filenames
wc -l $filenames

# Generate bad reads summary
cd ~/final_ngs/dnaseq/data/untrimmed_fastq
for filename in *.fastq.gz
do 
  echo $filename 
  zgrep -B1 -A2 NNNNNNNNNN $filename > $filename-bad-reads.fastq
  zgrep -cH NNNNNNNNNN $filename >> bad-reads.count.summary
done

# Move files
mkdir ../other
mv ~/final_ngs/dnaseq/data/untrimmed_fastq/*bad* ~/final_ngs/dnaseq/untrimmed_fastq/other

# Assessing quality
cd ~/final_ngs/dnaseq/data/untrimmed_fastq
fastqc -t 4 *.fastq.gz
mkdir ~/final_ngs/dnaseq/results/fastqc_untrimmed_reads
mv *fastqc* ~/final_ngs/dnaseq/results/fastqc_untrimmed_reads/
ls -lh ~/final_ngs/dnaseq/results/fastqc_untrimmed_reads/
# Get files in FileZilla and view  HTML

# Unzip other files
cd ~/final_ngs/dnaseq/results/fastqc_untrimmed_reads/
for zip in *.zip; do unzip $zip; done
ls -lh SOME_FILE_fastqc
head SOME_FILE_fastqc/summary.txt
cat */summary.txt > ~/final_ngs/dnaseq/logs/fastqc_summaries.txt

# Trimming
cd ~/final_ngs/dnaseq/data/untrimmed_fastq
java -jar /home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/trimmomatic.jar PE \
-threads 4 -phred33 \
/home/ubuntu/final_ngs/dnaseq/data/untrimmed_fastq/NGS0001.R1.fastq.gz \
/home/ubuntu/final_ngs/dnaseq/data/untrimmed_fastq/NGS0001.R2.fastq.gz \
-baseout /home/ubuntu/final_ngs/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R \
ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10 \
TRAILING:25 MINLEN:50

# FastQC on trimmed reads
fastqc -t 4 /home/ubuntu/final_ngs/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R_1P.fastq.gz \
/home/ubuntu/final_ngs/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R_1U.fastq.gz \
/home/ubuntu/final_ngs/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R_2P.fastq.gz \
/home/ubuntu/final_ngs/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R_2U.fastq.gz

# Fix file format
mv NGS0001_trimmed_R_1P NGS0001_trimmed_R_1P.fastq.gz
mv NGS0001_trimmed_R_1U NGS0001_trimmed_R_1U.fastq.gz
mv NGS0001_trimmed_R_2P NGS0001_trimmed_R_2P.fastq.gz
mv NGS0001_trimmed_R_2U NGS0001_trimmed_R_2U.fastq.gz

# BWA indexing
cd ~/final_ngs/dnaseq/data
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
mkdir -p ~/final_ngs/dnaseq/data/reference
mv ~/final_ngs/dnaseq/data/hg19.fa.gz ~/final_ngs/dnaseq/data/reference/
bwa index ~/final_ngs/dnaseq/data/reference/hg19.fa.gz

# Aligning
cd ~/final_ngs/dnaseq/data/trimmed_fastq
mv NGS0001_trimmed_R_1P.fastq.gz NGS0001_trimmed_R_1P.fastq
mv NGS0001_trimmed_R_2P.fastq.gz NGS0001_trimmed_R_2P.fastq
mkdir -p ~/final_ngs/dnaseq/data/aligned_data
bwa mem -t 4 -v 1 -R '@RG\tID:11V6WR1:111:D1375ACXX:1:1101:1671:2229\tSM:NGS0001\tPL:ILLUMINA\tLB:nextera-ngs0001\tDT:2024-04-05\tPU:HWI-D00119' \
-I 250,50 ~/final_ngs/dnaseq/data/reference/hg19.fa.gz \
~/final_ngs/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R_1P.fastq \
~/final_ngs/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R_2P.fastq \
> ~/final_ngs/dnaseq/data/aligned_data/NGS0001.sam

# Convert SAM to BAM
cd ~/final_ngs/dnaseq/data/aligned_data
samtools view -h -b NGS0001.sam > NGS0001.bam
samtools sort NGS0001.bam > NGS0001_sorted.bam
samtools index NGS0001_sorted.bam

# Mark duplicates
picard MarkDuplicates I=NGS0001_sorted.bam O=NGS0001_sorted_marked.bam M=marked_dup_metrics.txt
samtools index NGS0001_sorted_marked.bam

# Filter on mapping quality
samtools view -F 1796 -q 20 -o NGS0001_sorted_filtered.bam NGS0001_sorted_marked.bam
samtools index NGS0001_sorted_filtered.bam

# Flagstats
samtools flagstat NGS0001_sorted.bam > NGS0001_flagstats.txt
samtools idxstats NGS0001_sorted.bam > NGS0001_idxstats.txt
samtools depth NGS0001_sorted.bam > NGS0001_coverage.txt

#view summary results

head NGS0001_flagstats.txt
head  NGS0001_idxstats.txt
head NGS0001_coverage.txt


# Collect insert size metrics
java -jar /home/ubuntu/anaconda3/pkgs/picard-2.18.29-0/share/picard-2.18.29-0/picard.jar CollectInsertSizeMetrics \
  I=NGS0001_sorted.bam \
  O=NGS0001_insert_metrics.txt \
  H=NGS0001_insert_size_histogram.pdf

# Variant calling with FreeBayes
zcat ~/final_ngs/dnaseq/data/reference/hg19.fa.gz > ~/final_ngs/dnaseq/data/reference/hg19.fa
samtools faidx ~/final_ngs/dnaseq/data/reference/hg19.fa
freebayes --bam ~/final_ngs/dnaseq/data/aligned_data/NGS0001_sorted_filtered.bam \
  --fasta-reference ~/final_ngs/dnaseq/data/reference/hg19.fa \
  --vcf ~/final_ngs/dnaseq/results/NGS0001.vcf \
  -i "QUAL >= 20 && DP >= 10"
  
bgzip ~/final_ngs/dnaseq/results/NGS0001.vcf
tabix -p vcf ~/final_ngs/dnaseq/results/NGS0001.vcf.gz

# Annotation with ANNOVAR
cd ~/
tar -zxvf annovar.latest.tar.gz
cd annovar
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro humandb/
./convert2annovar.pl -format vcf4 ~/final_ngs/dnaseq/results/NGS0001.vcf.gz > ~/final_ngs/dnaseq/results/NGS0001.avinput
./table_annovar.pl ~/final_ngs/dnaseq/results/NGS0001.avinput humandb/ -buildver hg19 \
  -out ~/final_ngs/dnaseq/results/NGS0001 -remove \
  -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f \
  -otherinfo -nastring . -csvout
#download .csv from filezilla and view in excell

# Annotation with SnpEff
cd ~/final_ngs/dnaseq
wget http://downloads.sourceforge.net/project/snpeff/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
cd snpEff
java -jar snpEff.jar download -v hg19
java -Xmx4g -jar snpEff.jar eff -v -no-intergenic -i vcf -o vcf hg19 \
  ~/final_ngs/dnaseq/results/NGS0001.vcf.gz > ~/final_ngs/dnaseq/results/NGS0001_snpeff.vcf

# Open resulting HTML file with FileZilla


