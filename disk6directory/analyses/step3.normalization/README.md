# Step 3: Normalization of raw read counts

*A General Overview*: 

The purpose of normalizing counts is to remove *systematic variation* in sample read counts.  This variation can be due to fixed factors such as: time of sample collection, index used, or the physical sequencing machine.  For these samples, we have meta-data on index, lane, batch, effective transcript length, and pedigree relationships -- All of which can potentially affect or contribute to the magnitude of counts.

Given the objective is to enhance the detection of shared expression patterns between families who have similar phenotypes, it would nice to leverage a covariance matrix describing the relationships among samples.  This type of covariance matrix, along with the other systematic variation, can be taken into consideration when using a mixed effect model.  

Utilizing the `ASReml-R`  package we can model log2-transformed count data using a guassian distribution to remove the effects of index, lane, batch, and pedigree while attempting to estimate the random effect component for each biological replicate.

A mixed-model approach is different than modeling counts with a negative-binomial or poisson distribution, which is what is typically done with RNA-Seq reads.  Therefore as a comparison, we can utilize the R package `DESeq2` to normalize the biological replicate counts in the standard fashion.  Since read counts were generated with Salmon, the `tximport` package can be used in conjunction with DESeq to normalize reads while accounting for index, batch, and effective length estimates.


## Method and code

### Linear Mixed Model

A linear mixed-model was implemented with technical replicate count data using ASReml-R to normalize for: batch, lane, index, and pedigree.  

The model used was:


    asreml(fixed= txpt ~1 + batch + lane + index_seq, random = ~ ide(animal_id,var=T,init=1) + ped(fam_id,var=T,init=1),
    ginverse=list(fam_id=ainv,animal_id=ainv2),data=tip.subset,maxiter=50,trace=F, family=asreml.gaussian(link = "identity"))

To see the complete input and output head to: [LMM Normalization](https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/analyses/step3.normalization/LMM_Norm.Rmd)

Output is located at: `/media/disk6/ARF/RNASEQ/counts/86kSalmon/Step3_LMM_Counts.RData`

### DESeq2 normalization

DESeq2 was implemented with biological replicate count data to normalize for: batch, index, and effective length.

The model used was:

To see the complete input and output head to: 

Output is located at:
