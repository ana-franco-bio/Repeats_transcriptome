#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -l h_vmem=10G
#$ -l h_rt=6:0:0

module load samtools
module load bowtie2
#Chose clusters that have any percentage of mitochondria or plastid
#Go to RepeatExplorer output and open file: dna_database_annotation to see which reads came from organelles
#Go to reads.fasta and get them
#Select just interested reads
xargs samtools faidx reads35.fasta < reads_mit_35.txt > wanted_seq_only35.fasta
xargs samtools faidx reads75.fasta < reads_mit_75.txt > wanted_seq_only75.fasta
xargs samtools faidx reads72.fasta < reads_mit_72.txt > wanted_seq_only72.fasta
xargs samtools faidx reads60.fasta < reads_mit_60.txt > wanted_seq_only60.fasta
xargs samtools faidx reads56.fasta < reads_mit_56.txt > wanted_seq_only56.fasta
xargs samtools faidx reads97.fasta < reads_mit_97.txt > wanted_seq_only97.fasta

#Create index
bowtie2-build contigsCL35.fasta contigs35.ref
bowtie2-build contigsCL75.fasta contigs75.ref
bowtie2-build contigsCL72.fasta contigs72.ref
bowtie2-build contigsCL60.fasta contigs60.ref
bowtie2-build contigsCL56.fasta contigs56.ref
bowtie2-build contigsCL97.fasta contigs97.ref


#map reads to contigs
bowtie2 -f -x contigs75.ref wanted_seq_only35.fasta |samtools view -bS -> mapCL35_align.bam
bowtie2 -f -x contigs75.ref wanted_seq_only75.fasta |samtools view -bS -> mapCL75_align.bam
bowtie2 -f -x contigs72.ref wanted_seq_only72.fasta |samtools view -bS -> mapCL72_align.bam
bowtie2 -f -x contigs60.ref wanted_seq_only60.fasta |samtools view -bS -> mapCL60_align.bam
bowtie2 -f -x contigs56.ref wanted_seq_only56.fasta |samtools view -bS -> mapCL56_align.bam
bowtie2 -f -x contigs97.ref wanted_seq_only97.fasta |samtools view -bS -> mapCL97_align.bam

#Prepare bam files to igv
samtools sort mapCL35_align.bam -o mapCL35_align.bam.sorted.bam
samtools index mapCL35_align.bam.sorted.bam
samtools sort mapCL75_align.bam -o mapCL75_align.bam.sorted.bam
samtools index mapCL75_align.bam.sorted.bam
samtools sort mapCL72_align.bam -o mapCL72_align.bam.sorted.bam
samtools index mapCL72_align.bam.sorted.bam
samtools sort mapCL60_align.bam -o mapCL60_align.bam.sorted.bam
samtools index mapCL60_align.bam.sorted.bam
samtools sort mapCL56_align.bam -o mapCL56_align.bam.sorted.bam
samtools index mapCL56_align.bam.sorted.bam
samtools sort mapCL97_align.bam -o mapCL97_align.bam.sorted.bam
samtools index mapCL97_align.bam.sorted.bam
