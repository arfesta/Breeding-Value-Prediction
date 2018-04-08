# Breeding-Value-Prediction

*Alternative Methods to Breeding Value Prediction in Loblolly Pine*

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

Data | Path | Notes
--- | --- | ---
**raw read files** | `/media/disk6/ARF/RNASEQ/rawreads/86kSalmon` | Raw files returned from GSL
*raw tar* | `./EWtarfiles or ./LGEPtarfiles` | 
*raw fasta* | `./EWfasta or ./LGEPfasta` | 
**trimmed and filtered read files** | `/media/disk6/ARF/RNASEQ/trimmedfiltreads/86k` | Files post trim & adapater removal
*EW* | `./EW/lane01 ... ./lane12` | 
*LGEP* | `./LGEP/lane01 ... ./lane18` | 
**salmon count files** 
*EW tech reps* | `/media/disk6/ARF/RNASEQ/counts/86kSalmon/EW/` | `./lane01 ... ./lane12`
*LGEP tech reps* | `/media/disk6/ARF/RNASEQ/counts/86kSalmon/LGEP/` | `./lane01 ... ./lane18`
*EW bio reps* | `/media/disk6/ARF/RNASEQ/counts/86kSalmon/bio_EW/` | `./Sample_<animal_id>/`
*LGEP bio reps* | `/media/disk6/ARF/RNASEQ/counts/86kSalmon/bio_LGEP/` | `./Sample_<animal_id>/`
**experimental data resources**  | `/media/disk6/ARF/RNASEQ/Breeding-Value-Prediction/disk6directory/resources`
*sequencing* | `./exptdesign/sequencing` | `./EWtarfiles or ./LGEPtarfiles`
*pedigree* | `./pedigree` | `./EWfasta or ./LGEPfasta`
*phenotypes* | `./phenos` | `./EWfasta or ./LGEPfasta`
