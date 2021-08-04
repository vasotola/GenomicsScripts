#!/bin/bash
#SBATCH --job-name=parentFilter
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --mem=34gb
#SBATCH --time=167:00:00
#SBATCH --output=parentFilter.%j.out
#SBATCH --error=parentFilter.%j.err

#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=vas59108@uga.edu

###### IMPORTANT NOTE: ANYWHERE YOU SEE "cd" or "mv" or the reference genome YOU WILL NEED TO CHANGE THAT TO YOUR SPECIFIC DIRECTORIES

cd /scratch/vas59108/parCardRils_v2/testProtocol/parents/vcfFiltering/

module load GATK/4.1.6.0-GCCcore-8.3.0-Java-1.8

## change "DP < XX" to your min acceptable depth value
## change "DP > XX" to your mean+2*sd value of depth for this parent

## Parent 1 depth filter (par)
gatk --java-options "-Xmx4g" VariantFiltration  \
  -R /scratch/vas59108/references/cardRef/CE10g_v2.0.fa \
  -V YY1.12H.bam.g.vcf.gz \
  --filter-expression "DP < 5" \
  --filter-name "minDP_filter" \
  --filter-expression "DP > 21" \
  --filter-name "maxDP_filter" \
  -O par.DPfilt.vcf.gz 

## Parent 2 depth filter (card)
gatk --java-options "-Xmx4g" VariantFiltration  \
  -R /scratch/vas59108/references/cardRef/CE10g_v2.0.fa \
  -V YY1.01A.bam.g.vcf.gz \
  --filter-expression "DP < 5" \
  --filter-name "minDP_filter" \
  --filter-expression "DP > 25" \
  --filter-name "maxDP_filter" \
  -O card.DPfilt.vcf.gz 

## Combine files
gatk --java-options "-Xmx4g"  CombineGVCFs \
   -R /scratch/vas59108/references/cardRef/CE10g_v2.0.fa \
    --variant par.DPfilt.vcf.gz     \
    --variant card.DPfilt.vcf.gz    \
   -O parentsDP.vcf.gz

## Joint genotype 
gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R /scratch/vas59108/references/cardRef/CE10g_v2.0.fa \
   -V parentsDP.vcf.gz \
   -O gtype_parents.vcf.gz

## Quality filters
gatk --java-options "-Xmx4g" VariantFiltration  \
  -R /scratch/vas59108/references/cardRef/CE10g_v2.0.fa \
  -V gtype_parents.vcf.gz \
--filter-expression "QD < 2.0" \
--filter-name "QD_filter" \
--filter-expression "FS > 60.0" \
--filter-name "FS_filter" \
--filter-expression "MQ < 40.0" \
--filter-name "MQ_filter" \
--filter-expression "MQRankSum < -12.5" \
--filter-name "MQRankSum_filter" \
--filter-expression "ReadPosRankSum < -8.0" \
--filter-name "ReadPosRankSum_filter" \
  -O gtype_parents_qualDPfilt.vcf.gz
    
## Select SNPs and biallelic sites
gatk --java-options "-Xmx4g" SelectVariants  \
  -R /scratch/vas59108/references/cardRef/CE10g_v2.0.fa \
  -V gtype_parents_qualDPfilt.vcf.gz \
  --select-type-to-include SNP \
  --restrict-alleles-to BIALLELIC \
  -O gtype_parents_SNPs.vcf.gz
