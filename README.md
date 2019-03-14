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
    + [Step 6 - Prediction using LOO](#step-6---prediction-using-loo)

## Abstract

* Phenotypic variation in forest trees may be partitioned into genomic and environmntal compenets which are consequently used to estimate the heritability of traits as the proportion of total phenotypic variation attributed to genetic variation.

* Applied tree breeding programs can use matrices of relationships, based either on recorded pedigrees in structured breeding populations or on genotypes of molecular genetic markers, to model genetic covariation among related individuals and predict genetic values for individuals for whom no phenotypic measurements are available. 

* ***This study tests the hypothesis that genetic covariation among individuals of similar genetic value will be reflected in shared patterns of gene expression.*** We collected gene expression data by high-throughput sequencing of RNA isolated from pooled seedlings of parents with known genetic value, and compared alternative approaches to data analysis to test this hypothesis.

## Background of samples

* All information about samples is located in the RNA Seq Data repo

### Step 1 - Transcript normalization and SNP filtering

The counts were normarlized *multiple* ways, however the following way was used for prediction:

  1.  Using **techincal replicate** counts and asreml to normalize for batch, index, lane, and pedigree:
  
  [markdown](http://htmlpreview.github.com/?https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/analyses/step3.normalization/Step3_LMM_animal.html)
  
  2.  For more normalization schemes using DESEQ2, edgeR, sommer in bio and tech see repo folder step3.normalization
  
  3. SNP's were filtered multiple ways: 
  
  [markdown](http://htmlpreview.github.com/?https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/analyses/step3.normalization/snp_filtering.html)
  
  4. The final data sets used for prediction were restructured to be in identical order: 
  
  [markdown](http://htmlpreview.github.com/?https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/analyses/step4.prediction/step1_data_restructure.html)
  

### Step 2 - Prediction of EW families with LGEP

* The EW vs. LGEP

#### Organizing test and train data sets

* EW and LGEP families were subset into train and test objects

  [markdown](http://htmlpreview.github.com/?https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/analyses/step4.prediction/step2_lgep.vs.ew_data.html)

#### Conduct prediction on EW

   Family mean estimates of counts and snps were used for prediction with OmicKriging and glmnet (lasso/ridge): 
   
   [EW predictions](http://htmlpreview.github.com/?https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/analyses/step4.prediction/step3_lgep.vs.ew_pred.v3.html)


### Step 5 - Prediction of 70-fold CV

#### Construct 70 test groups

 Instead of predicting across batch, here we split the complete data set into a 7-fold CV (repeated 10 times).  The cv groups were split so that each test fold had individuals which were spread across the phenotypic range: 
 
 [create 70 fold](http://htmlpreview.github.com/?https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/analyses/step4.prediction/step4_create_70fold.html)

#### Conduct prediction on each of the test folds using all data

[prediction_script](https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/analyses/step4.prediction/step5_cv.70.predictions.V2.R)

#### Visualize predictions

[70 Fold CV visualization](https://htmlpreview.github.com/?https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/analyses/step4.prediction/step6_visualize_predictions.v2.html)


### Step 6 - Prediction using LOO

 Predictions were conducted using a maximum training size of 55 to predict the 56th family using OK, lasso, & ridge. Script: 
 
 [LOO Script](https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/analyses/step4.prediction/step7_loo.predictions.R)
  
#### Visualize prediction of LOO

[LOO markdown](https://htmlpreview.github.com/?https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/analyses/step4.prediction/step7_visualize_predictions.html)


## 
**The below part is defunct, the scripts are still there but are not used. 

 Estimate anova scores for features using LGEP and then conduct prediction on EW

 #### Generate anova scores using LGEP as training

  Utilizing the biological replicate data sets, ANOVA scores were estimated for each feature (snp/transcript): 
  
  [LGEP_ANOVA](http://htmlpreview.github.com/?https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/analyses/step4.prediction/step2_lgep.vs.ew_anova.html)
  
 #### Conduct prediction on EW

   Family mean estimates of counts and snps were used for prediction with OmicKriging and glmnet (lasso/ridge): 
   
   [EW predictions](http://htmlpreview.github.com/?https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/analyses/step4.prediction/step3_lgep.vs.ew_pred.html)
   
  

 **The below part is defunct, the scripts are still there but are not used.** 
 
  First construct the 70 test groups, estimate ANOVA scores, and then conduct predictions
 
 #### Conduct prediction on each of the test folds across pvals

Just as when predicting on the EW families, predictions were carried out for each of the 70 unique test groups: 

[predict 70-fold](https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/analyses/step4.prediction/step5_cv.70.predictions.Rmd)
  
 #### Visualize prediction of 70-fold

[70-fold-cv markdown](https://htmlpreview.github.com/?https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/analyses/step4.prediction/step6_visualize_predictions.html)

 
