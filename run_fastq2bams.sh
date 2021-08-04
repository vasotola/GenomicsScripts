#!/bin/bash
#SBATCH --job-name=fastq2bam
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --mem=32gb
#SBATCH --time=167:00:00
#SBATCH --output=fastq2bam.%j.out
#SBATCH --error=fastq2bam.%j.err

#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=vas59108@uga.edu

###### IMPORTANT NOTE: ANYWHERE YOU SEE "cd" or "mv" or the reference genome YOU WILL NEED TO CHANGE THAT TO YOUR SPECIFIC DIRECTORIES


module load BWA/0.7.17-GCC-8.3.0
module load Trimmomatic/0.39-Java-1.8.0_144
module load SAMtools/1.10-iccifort-2019.5.281


### Step 1: Reference Index for alignment

cd /scratch/vas59108/references/

bwa index CE10g_v2.0.fa
samtools faidx CE10g_v2.0.fa

### Step 2: unzip fastq files (if zipped)

cd /scratch/vas59108/parCardRils_v2/testProtocol/fastqFiles

gunzip *.gz

### Step 3: Remove adapters and low quality sequences

for i in *1.fq; do name=$(echo $i | rev | cut -c6- | rev)
	a=$name".1.fq"
	b=$name".2.fq"

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
		-threads 4 "${a}" "${b}" \
		"${a}.1.paired.fastq" \
		"${a}.1.unpaired.fastq" \
		"${b}.2.paired.fastq" \
		"${b}.2.unpaired.fastq" \
	     ILLUMINACLIP:paired_end.fa:2:20:10:4 \
	     LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

#This will perform the following:
	#Remove adapters (ILLUMINACLIP:TruSeq3-PE.fa:2:30:10)
	#Remove leading low quality or N bases (below quality 3) (LEADING:3)
	#Remove trailing low quality or N bases (below quality 3) (TRAILING:3)
	#Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15 (SLIDINGWINDOW:4:15)
	#Drop reads below the 36 bases long (MINLEN:36)

mv *unpaired.fastq /scratch/vas59108/parCardRils_v2/testProtocol/fastqFiles/unpairedData/
mv *paired.fastq /scratch/vas59108/parCardRils_v2/testProtocol/fastqFiles/pairedData/
#mv *.1.fq *.2.fq /scratch/vas59108/parCardRils_v2/testProtocol/fastqFiles/ogFastqs/

### Step 4: Map reads to genome (bwa)

### change directory of reference genome ###

cd /scratch/vas59108/parCardRils_v2/testProtocol/fastqFiles/pairedData/

for i in `ls *.1.fq.1.paired.fastq`; do newname=$(echo $i | rev | cut -c21- | rev)
	a=$newname".1.fq.1.paired.fastq"
	b=$newname".2.fq.2.paired.fastq"
	bwa mem /scratch/vas59108/references/cardRef/CE10g_v2.0.fa -t 8 "${a}" "${b}" > "${newname}.sam"
done

mv *.sam /scratch/vas59108/parCardRils_v2/testProtocol/samFiles/

### Step 5: Quality filter and sort sam, convert to bam

cd /scratch/vas59108/parCardRils_v2/testProtocol/samFiles

for f in *.sam
do
	samtools view -q 29 -b $f > ${f/.sam/}.bam
done

for f in *.bam
do
        samtools sort $f -o ${f%.*}_sorted.bam
done

mv *.bam /scratch/vas59108/parCardRils_v2/testProtocol/bamFiles

### Step 6: Add read groups

cd /scratch/vas59108/parCardRils_v2/testProtocol/bamFiles/

module load picard/2.16.0-Java-1.8.0_144

for f in *_sorted.bam
do
time java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
                I=$f \
                O=${f/_sorted.bam/}_RGsorted.bam \
                RGSM=$f \
                RGLB=NEBNext \
                RGPL=Illumina \
                RGPU=UNKNOWN VALIDATION_STRINGENCY=LENIENT; # adds read groups
done

### Step 7: Index bam files

cd /scratch/vas59108/parCardRils_v2/testProtocol/bamFiles/

for f in *RGsorted.bam
do 
	samtools index $f;
done
