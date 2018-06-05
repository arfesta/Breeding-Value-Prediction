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

```
After filtering, kept 192 out of 192 Individuals
Outputting VCF file...
After filtering, kept 663796 out of a possible 2269592 Sites
Run Time = 545.00 seconds
```

This filter kept 660K SNPs out of 2.2 million.

## Filtering minimum depth

The next filter we will apply is a minimum depth for a genotype call and a minimum mean depth. This command will recode genotypes that have less than 3 reads.

```
vcftools --vcf raw.g5mac3.recode.vcf --minDP 3 --recode --recode-INFO-all --out raw.g5mac3dp3
```

This step does not remove any SNPs

```
After filtering, kept 192 out of 192 Individuals
Outputting VCF file...
After filtering, kept 663796 out of a possible 663796 Sites
Run Time = 469.00 seconds
```

## Download error evaluation script

This script was found on the tutorial and can be used to estimate potential errors due to low read depth.

```
curl -L -O https://github.com/jpuritz/dDocent/raw/master/scripts/ErrorCount.sh
chmod +x ErrorCount.sh 
./ErrorCount.sh raw.g5mac3dp3.recode.vcf
```

Output:
```
This script counts the number of potential genotyping errors due to low read depth
It report a low range, based on a 50% binomial probability of observing the second allele in a heterozygote and a high range based on a 25% probability.
```

## Assess the level of missing data

```
vcftools --vcf raw.g5mac3dp3.recode.vcf --missing-indv
```

This will create an output called out.imiss. Let’s examine it.

```
cat out.imiss
```

You can see that some individuals have as high as 61% missing data. We definitely want to filter those out. Let’s take a look at a histogram

```
mawk '!/IN/' out.imiss | cut -f5 > totalmissing
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of % missing data per individual"
set ylabel "Number of Occurrences"
set xlabel "% of missing data"
#set yr [0:100000]
binwidth=0.01
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'totalmissing' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF
```

Now make a list of the individuals which are missing 50% of the data

```
mawk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv
```
INDV
85.ew.bio.rep.merge.bam
9_.ew.bio.rep.merge.bam

Remove these two individuals

```
vcftools --vcf raw.g5mac3dp3.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out raw.g5mac3dplm
```
