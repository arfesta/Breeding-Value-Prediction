# Breeding-Value-Prediction

*Alternative Methods to Breeding Value Prediction in Loblolly Pine*

- [Breeding-Value-Prediction](#breeding-value-prediction)
  * [Abstract](#abstract)
  * [Background of samples](#background-of-samples)
  * [Location of data](#location-of-data)
  * [Analyses](#analyses)
    + [Step 1 - Data Prep](#step-1---data-prep)
    + [Step 2 - Load Count Data](#step-2---load-count-data)
    + [Step 3 - Transcript normalization and SNP filtering](#step-3---transcript-normalization-and-snp-filtering)
    + [Step 4 - Prediction of EW families with LGEP](#step-4---prediction-of-ew-families-with-lgep)
    + [Step 5 - Prediction of 70-fold CV](#step-5---prediction-of-70-fold-cv)
    
## Abstract

* Phenotypic variation in forest trees may be partitioned into genomic and environmntal compenets which are consequently used to estimate the heritability of traits as the proportion of total phenotypic variation attributed to genetic variation.

* Applied tree breeding programs can use matrices of relationships, based either on recorded pedigrees in structured breeding populations or on genotypes of molecular genetic markers, to model genetic covariation among related individuals and predict genetic values for individuals for whom no phenotypic measurements are available. 

* ***This study tests the hypothesis that genetic covariation among individuals of similar genetic value will be reflected in shared patterns of gene expression.*** We collected gene expression data by high-throughput sequencing of RNA isolated from pooled seedlings of parents with known genetic value, and compared alternative approaches to data analysis to test this hypothesis.

## Background of samples

* LGEP
   
   - 144 Biological Replicates x 3 technical replicates = 432 technical replicates

* EW

   - 80 Biological Replicates x 3 technical replicates = 240 technical replicates

## Location of data

Data Subject Type | Data File Type | Path | Notes
--- | --- | --- | ---
**raw read files**  | | `/media/disk6/ARF/RNASEQ/shared/rawreads/86kSalmon` | Raw files returned from GSL
| | *raw tar* | `./EWtarfiles or ./LGEPtarfiles`
| | *raw fasta* | `./EWfasta or ./LGEPfasta` 
**trimmed and filtered read files** | |`/media/disk6/ARF/RNASEQ/shared/trimmedfiltreads/86k` | Files post trim & adapater removal
|  |*EW* | `./EW/lane01 ... ./lane12` | 
|  |*LGEP* | `./LGEP/lane01 ... ./lane18` | 
**salmon count files** | |`/media/disk6/ARF/RNASEQ/shared/counts/86kSalmon` | Direcotries containing quant.sf files
|  |*EW tech reps* | `./EW/lane01 ... ./lane12` | 
|  |*LGEP tech reps* | `./LGEP/lane01 ... ./lane18` |
|  |*EW bio reps* | `./bio_EW/Sample_<animal_id>/` | 
|  |*LGEP bio reps* | `./bio_LGEP/Sample_<animal_id>/` | 
**experimental data resources**  | | `/media/disk6/ARF/RNASEQ/BV-Prediction/Breeding-Value-Prediction/disk6directory/resources` | Experiment information
|  |*sequencing* | `./exptdesign/sequencing` | `./EWtarfiles or ./LGEPtarfiles`
|  |*pedigree* | `./pedigree` | `./EWfasta or ./LGEPfasta`
|  |*phenotypes* | `./phenos` | `./EWfasta or ./LGEPfasta`


## Analyses

### Step 1 - Data Prep

   Data prep includes everything from unpacking the original tar files recieved by the GSL up to estimating transcript               
      abundance with Salmon. Additionally, this step includes identification of the indicies used within both batches and creates an experimental info matrix containing all meta data from both batches.
      
   See the [raw reads README](http://htmlpreview.github.com/?https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/rawreads/012718raw_data_processing.html) for step by step processing of files.  
   
  Extra prep scripts: [sample index identification](http://htmlpreview.github.com/?https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/analyses/step1.data_prep/identify_index_used.ouput.html) & [creation of experimental data](http://htmlpreview.github.com/?https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/analyses/step1.data_prep/create_expt_data.html)

### Step 2 - Load Count Data
      
   Once counts have been estimated, the next step involves reading in the aligned technical, or biological, replicate counts using the tximport package.
      
   Additionally, the phenotype and other sample meta-data is constructed for normalization.
   
   To see this process for **biological reps**, navigate to: [load counts bio rep html file](http://htmlpreview.github.com/?https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/analyses/step2.loadcounts/load.counts.html) which contains the complete markdown and output.

   To see this process for **technical reps**, navigate to: [load counts tech rep html file](http://htmlpreview.github.com/?https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/analyses/step2.loadcounts/load.counts_techreps.html) which contains the complete markdown and output.

### Step 3 - Transcript normalization and SNP filtering

The counts were normarlized *multiple* ways, however only the #1 was used for prediction:

  1.  Using **techincal replicate** counts and asreml to normalize for batch, index, lane, and pedigree: [markdown](http://htmlpreview.github.com/?https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/analyses/step3.normalization/Step3_LMM_animal.html)
  
  2.  Using DESEQ2, edgeR, sommer in bio and tech see repo
  
  3. SNP's were filtered multiple ways: [markdown](http://htmlpreview.github.com/?https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/analyses/step3.normalization/snp_filtering.html)
  
  4. The final data sets used for prediction were restructured to be in identical order: [markdown](http://htmlpreview.github.com/?https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/analyses/step4.prediction/step1_data_restructure.html)
  

### Step 4 - Prediction of EW families with LGEP

#### Conduct prediction on EW

   Family mean estimates of counts and snps were used for prediction with OmicKriging and glmnet (lasso/ridge): [EW predictions](http://htmlpreview.github.com/?https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/analyses/step4.prediction/step3_lgep.vs.ew_pred.v2.html)


#### **The below part is defunct, the scripts are still there but are not used. 

 Estimate anova scores for features using LGEP and then conduct prediction on EW

#### Generate anova scores using LGEP as training

  Utilizing the biological replicate data sets, ANOVA scores were estimated for each feature (snp/transcript): [LGEP_ANOVA](http://htmlpreview.github.com/?https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/analyses/step4.prediction/step2_lgep.vs.ew_anova.html)
  
#### Conduct prediction on EW

   Family mean estimates of counts and snps were used for prediction with OmicKriging and glmnet (lasso/ridge): [EW predictions](http://htmlpreview.github.com/?https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/analyses/step4.prediction/step3_lgep.vs.ew_pred.html)
   
  
### Step 5 - Prediction of 70-fold CV

 First construct the 70 test groups, estimate ANOVA scores, and then conduct predictions

#### Construct 70 test groups

 Instead of predicting across batch, here we split the complete data set into a 7-fold CV (repeated 10 times).  The cv groups were split so that each test fold had individuals which were spread across the phenotypic range: [create 70 fold](http://htmlpreview.github.com/?https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/analyses/step4.prediction/step4_create_70fold.html)
 
#### Conduct prediction on each of the test folds

Just as when predicting on the EW families, predictions were carried out for each of the 70 unique test groups: [predict 70-fold](https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/analyses/step4.prediction/step5_cv.70.predictions.Rmd)
  
#### Visualize prediction of 70-fold

[70-fold-cv markdown](https://htmlpreview.github.com/?https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/analyses/step4.prediction/step6_visualize_predictions.html)
 
 
