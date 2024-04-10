#!/bin/bash

# M2305834 NGS pipeline 

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

# Indexing reference genome
cd ~/final_ngs/dnaseq/data
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
mkdir -p ~/final_ngs/dnaseq/data/reference
mv ~/final_ngs/dnaseq/data/hg19.fa.gz ~/final_ngs/dnaseq/data/reference/
bwa index ~/final_ngs/dnaseq/data/reference/hg19.fa.gz

# Variant calling with GATK
cd ~/final_ngs/dnaseq/data/aligned_data

java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R ~/final_ngs/dnaseq/data/reference/hg19.fa \
                                 -I NGS0001_sorted_filtered.bam \
                                 -o NGS0001_gatk.vcf

bgzip NGS0001_gatk.vcf
tabix -p vcf NGS0001_gatk.vcf.gz

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
./convert2annovar.pl -format vcf4 ~/final_ngs/dnaseq/results/NGS0001_gatk.vcf.gz > ~/final_ngs/dnaseq/results/NGS0001_gatk.avinput
./table_annovar.pl ~/final_ngs/dnaseq/results/NGS0001_gatk.avinput humandb/ -buildver hg19 \
  -out ~/final_ngs/dnaseq/results/NGS0001_gatk -remove \
  -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f \
  -otherinfo -nastring . -csvout

# Annotation with SnpEff
cd ~/final_ngs/dnaseq
wget http://downloads.sourceforge.net/project/snpeff/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
cd snpEff
java -jar snpEff.jar download -v hg19
java -Xmx4g -jar snpEff.jar eff -v -no-intergenic -i vcf -o vcf hg19 \
  ~/final_ngs/dnaseq/results/NGS0001_gatk.vcf.gz > ~/final_ngs/dnaseq/results/NGS0001_gatk_snpeff.vcf

# Open resulting HTML file with FileZilla

