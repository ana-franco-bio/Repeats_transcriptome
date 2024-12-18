#!/bin/sh
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -l h_rt=12:0:0
#$ -l h_vmem=20G

###input reads

INPUT_R1=Setaria2-ribominus_S28_R1_001.fastq.gz
INPUT_R2=Setaria2-ribominus_S28_R2_001.fastq.gz
TRANSLIB=Setaria2-ribominus

MINlENGTH=150 #minimus length for trimmomatic length filtering

module load java
#trimmomatic length filtering
java -jar /data/home/btx671/programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE $INPUT_R1 $INPUT_R2 trimmed${INPUT_R1}_paired.fastq trimmed${INPUT_R1}_unpaired.fastq trimmed${INPUT_R2}_paired.fastq trimmed${INPUT_R2}_unpaired.fastq ILLUMINACLIP:/data/home/btx671/programs/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:$MINlENGTH HEADCROP:0

#run Fastqc again
module unload java
module load fastqc

fastqc trimmed${INPUT_R1}_paired.fastq
fastqc trimmed${INPUT_R2}_paired.fastq


#load the modules
module load samtools
module load bowtie2

#create reference
bowtie2-build allrepeats.fasta allrepeats.ref

#bowtie2 map
bowtie2 -q -x allrepeats.ref -1 trimmed${INPUT_R1}_paired.fastq -2 trimmed${INPUT_R2}_paired.fastq|samtools view -bS -> lengthfiltered_allrepeats_${TRANSLIB}_align.bam

#samtools filter out the unmapped reads
samtools view -h -F 4 -b lengthfiltered_allrepeats_${TRANSLIB}_align.bam > mappedlengthfiltered_allrepeats_${TRANSLIB}_align.bam
samtools sort mappedlengthfiltered_allrepeats_${TRANSLIB}_align.bam -o mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted.bam
samtools view mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted.bam|cut -f1,3 > mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt

#Loop the reads number for each cluster
for i in {1..208}
do
  grep "CL${i}Contig" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt | wc -l >> file${TRANSLIB}.txt
 done
 
#organeles 
grep -E "CL3Contig8" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}file_Reads_Mitochondria_CL3.txt
grep -E "CL3Contig38" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}file_Reads_Mitochondria_CL3.txt
grep -E "CL3Contig47" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}file_Reads_Mitochondria_CL3.txt
grep -E "CL3Contig63" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}file_Reads_Mitochondria_CL3.txt
grep -E "CL3Contig140" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}file_Reads_Mitochondria_CL3.txt
grep -E "CL3Contig205" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}file_Reads_Mitochondria_CL3.txt
grep -E "CL3Contig210" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}file_Reads_Mitochondria_CL3.txt
grep -E "CL3Contig224" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}file_Reads_Mitochondria_CL3.txt
grep -E "CL3Contig311" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}file_Reads_Mitochondria_CL3.txt
grep -E "CL3Contig312" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}file_Reads_Mitochondria_CL3.txt
grep -E "CL3Contig319" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}file_Reads_Mitochondria_CL3.txt
grep -E "CL3Contig352" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}file_Reads_Mitochondria_CL3.txt
grep -E "CL3Contig401" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}file_Reads_Mitochondria_CL3.txt

#CL35
grep -E "CL35Contig78" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}file_Reads_Mitochondria_CL35.txt
grep -E "CL35Contig115" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}file_Reads_Mitochondria_CL35.txt
grep -E "CL35Contig43" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}file_Reads_Mitochondria_CL35.txt
grep -E "CL35Contig76" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}file_Reads_Mitochondria_CL35.txt

#CL56
grep -E "CL56Contig26" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}file_Reads_Mitochondria_CL56.txt
grep -E "CL56Contig33" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}file_Reads_Mitochondria_CL56.txt

#CL60
grep -E "CL60Contig45" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}file_Reads_Mitochondria_CL60.txt

#CL75
grep -E "CL75Contig11" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}file_Reads_Mitochondria_CL75.txt

#CL95
grep -E "CL95Contig5" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}file_Reads_Mitochondria_CL95.txt

#CL72
grep -E "CL72Contig72" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}file_Reads_Mitochondria_CL72.txt
grep -E "CL72Contig124" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}file_Reads_Mitochondria_CL124.txt

#Clusters with tRNA

#CL2
grep -E "CL2Contig40" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL2.txt
grep -E "CL2Contig59" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL2.txt
grep -E "CL2Contig69" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL2.txt
grep -E "CL2Contig495" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL2.txt
grep -E "CL2Contig612" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL2.txt
grep -E "CL2Contig731" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL2.txt
grep -E "CL2Contig735" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL2.txt
grep -E "CL2Contig752" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL2.txt


#CL7
grep -E "CL7Contig16" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL7.txt
grep -E "CL7Contig24" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL7.txt
grep -E "CL7Contig78" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL7.txt
grep -E "CL7Contig99" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL7.txt
grep -E "CL7Contig138" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL7.txt
grep -E "CL7Contig140" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL7.txt
grep -E "CL7Contig207" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL7.txt
grep -E "CL7Contig233" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL7.txt
grep -E "CL7Contig318" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL7.txt
grep -E "CL7Contig324" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL7.txt
grep -E "CL7Contig332" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL7.txt
grep -E "CL7Contig345" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL7.txt

#CL8
grep -E "CL8Contig67" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL8.txt
grep -E "CL8Contig145" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL8.txt
grep -E "CL8Contig160" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL8.txt
grep -E "CL8Contig194" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL8.txt
grep -E "CL8Contig203" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL8.txt
grep -E "CL8Contig239" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL8.txt
grep -E "CL8Contig306" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL8.txt
grep -E "CL8Contig316" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL8.txt
grep -E "CL8Contig362" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL8.txt
grep -E "CL8Contig480" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL8.txt
grep -E "CL8Contig523" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL8.txt
grep -E "CL8Contig613" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL8.txt
grep -E "CL8Contig735" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL8.txt
grep -E "CL8Contig794" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL8.txt

#CL18
grep -E "CL18Contig92" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL18.txt
grep -E "CL18Contig184" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL18.txt
grep -E "CL18Contig305" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL18.txt

#CL19
grep -E "CL19Contig160" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL19.txt
grep -E "CL19Contig214" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL19.txt
grep -E "CL19Contig277" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL19.txt
grep -E "CL19Contig444" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL19.txt

#CL22
grep -E "CL22Contig207" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL22.txt
grep -E "CL22Contig239" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL22.txt
grep -E "CL22Contig351" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL22.txt

#CL25
grep -E "CL25Contig131" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL25.txt
grep -E "CL25Contig192" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL25.txt
grep -E "CL25Contig298" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL25.txt
grep -E "CL25Contig318" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL25.txt

#CL26
grep -E "CL26Contig92" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL26.txt
grep -E "CL26Contig96" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL26.txt

#CL40
grep -E "CL40Contig16" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL40.txt
grep -E "CL40Contig59" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL40.txt
grep -E "CL40Contig62" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL40.txt
grep -E "CL40Contig68" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL40.txt
grep -E "CL40Contig72" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL40.txt
grep -E "CL40Contig98" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL40.txt

#CL41
grep -E "CL41Contig300" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL41.txt

#CL60
grep -E "CL60Contig321" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL60.txt

#CL81
grep -E "CL81Contig36" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL81.txt
grep -E "CL81Contig45" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL81.txt

#CL99
grep -E "CL99Contig7" mappedlengthfiltered_allrepeats_${TRANSLIB}_align_sorted_cut.txt| wc -l >> ${TRANSLIB}_tRNA_CL99.txt