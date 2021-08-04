#!/bin/bash
#SBATCH --job-name=bam2gvcf
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --mem=32gb
#SBATCH --time=167:00:00
#SBATCH --output=bam2gvcf.%j.out
#SBATCH --error=bam2gvcf.%j.err

#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=vas59108@uga.edu

###### IMPORTANT NOTE: ANYWHERE YOU SEE "cd" or "mv" or the reference genome YOU WILL NEED TO CHANGE THAT TO YOUR SPECIFIC DIRECTORIES

module load GATK/4.1.6.0-GCCcore-8.3.0-Java-1.8

### Step 1: Make a .dict file (required by gatk)

cd /scratch/vas59108/references/cardRef/

gatk CreateSequenceDictionary -R CE10g_v2.0.fa

### Step 2: convert bams to gvcf files

cd /scratch/vas59108/parCardRils_v2/testProtocol/bamFiles/

for i in *.bam;
do
gatk --java-options "-Xmx4g" HaplotypeCaller  \
  -R /scratch/vas59108/references/cardRef/CE10g_v2.0.fa \
  -I $i \
  -O $i.g.vcf.gz \
  -ERC GVCF
done
