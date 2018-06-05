# SNP Filtering

Starting at the end of the previous [SNP processing file](https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/snps/012918processing_snps_notes.Rmd), the purpose of this document is to practice filtering the 192 bioloigcal replicate merged SNP data for the purpose of prediction.

## Document set-up

```
cd /media/disk6/ARF/shared/snps
bcftools index ./192.biorep.vcf.bgzip.gz
mkdir filtering
cd filtering
```

## Filtering step 1:

The raw vcf gzipped file (192.biorep.vcf.bgzip.gz) is going to have a lot of erroneous variant calls and a lot of variants that are only present in one individual.

To make this file more manageable, let’s start by applying three step filter. We are going to only keep variants that have been successfully genotyped in 50% of individuals, a minimum quality score of 30, and a minor allele count of 3.

```
vcftools --gzvcf ../192.biorep.vcf.bgzip.gz --max-missing 0.5 --mac 3 --minQ 30 --recode --recode-INFO-all --out raw.g5mac3
```

In this code, we call vcftools, feed it a vcf file after the --gzvcf flag, --max-missing 0.5 tells it to filter genotypes called below 50% (across all individuals) the --mac 3 flag tells it to filter SNPs that have a minor allele count less than 3.
This is relative to genotypes, so it has to be called in at least 1 homozygote and 1 heterozygote or 3 heterozygotes. The --recode flag tells the program to write a new vcf file with the filters, --recode-INFO-all keeps all the INFO flags from the old vcf file in the new one.

## Filtering minimum depth

The next filter we will apply is a minimum depth for a genotype call and a minimum mean depth. This command will recode genotypes that have less than 3 reads.

```
vcftools --vcf raw.g5mac3.recode.vcf --minDP 3 --recode --recode-INFO-all --out raw.g5mac3dp3
```

## Download error evaluation script

```
curl -L -O https://github.com/jpuritz/dDocent/raw/master/scripts/ErrorCount.sh
chmod +x ErrorCount.sh 
./ErrorCount.sh raw.g5mac3dp3.recode.vcf
```
