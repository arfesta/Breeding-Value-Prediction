-   [4/13/18 Creating counts from raw
    reads](#creating-counts-from-raw-reads)
    -   [1. Extract raw data](#extract-raw-data)
        -   [Location of raw data](#location-of-raw-data)
        -   [Extracting LGEP tar files](#extracting-lgep-tar-files)
        -   [Extracting EW tar files](#extracting-ew-tar-files)
        -   [Untar sequencing lanes](#untar-sequencing-lanes)
    -   [2. Run fastqc on raw files](#run-fastqc-on-raw-files)
    -   [3. Filter and trim raw reads](#filter-and-trim-raw-reads)
    -   [4. FastQC can be re-done on filtered and trimmed
        files](#fastqc-can-be-re-done-on-filtered-and-trimmed-files)
    -   [5. Create complete experimental design
        matrix](#create-complete-experimental-design-matrix)
        -   [Check sample indicies](#check-sample-indicies)
        -   [Integrate phenotypes, pedigree, and sequencing
            data](#integrate-phenotypes-pedigree-and-sequencing-data)
    -   [6. Align technical replicates with
        Salmon](#align-technical-replicates-with-salmon)
    -   [7. Align biological replicates with
        Salmon](#align-biological-replicates-with-salmon)

4/13/18 Creating counts from raw reads
======================================

1. Extract raw data
-------------------

### Location of raw data

-   Raw mRNA sequencing libraries returned from the GSL are located at:

        /media/disk6/ARF/RNASEQ/rawreads/EWtarfiles && LGEPtarfiles

### Extracting LGEP tar files

\*Within the LGEP directory there are 18 .tar files (144 biological
samples) and then an additional directory which includes tar files from
the first run of sequencing. The 18 tar files were delivered at two
seperate times:

\*First set included 6 tar files labeled:

      -Project_20150504-GSL028B-Whetten-X1POOL && X2,X3,X4,X5,X6
      

\*Second set included 12 tar files labeled:

      -Project_20150615-GSL033A-Whetten-Y1POOL && 2,3,4,5,6 && Z1,Z2,Z3,Z4,Z5,Z6

\*The names of directores where manually edited to be:

      -lane01:lane06 corresponding to X1:X6 
      
      -lane07:lane12 corresponding to Y1:Y6
      
      -lane13:lane18 corresponding to Z1:Z6
     

\*Note sample in lane12 tar file 11x really belongs to lane11 so it was
moved there.

### Extracting EW tar files

\*Within the EW directory there are 12 .tar files (96 biological
samples). The 12 tar files were delivered at two seperate times:

\*First set included 4 tar files labeled:

      -20151203-GSL-046-B-Whetten-Pool9 && 10,11,12
     

\*Second set included 8 tar files labeled:

      -GSL-047-A-Whetten-Pool1 && 2,3,4,5,6,7,8

\*The names of directories were manually edited to be

      -lane01:lane12 corresponding to Pool  

### Untar sequencing lanes

-   List name of tar files and untar to specified directory

        cd /media/disk6/ARF/RNASEQ/rawreads/EWtarfiles
        ls *.tar > EW.tar.list

        while read file; do tar -xf ./$file --directory ../EWfasta/; done < EW.tar.list

        cd /media/disk6/ARF/RNASEQ/rawreads/LGEPtarfiles
        ls *.tar > LGEP.tar.list

        while read file; do tar -xf ./$file --directory ../LGEPfasta/; done < LGEP.tar.list

-   lanes 11,12,13,14,15,16 have folders of 24 samples with 3 smaller
    .fasta.gz files within each folder corresponding to that sample.
    These samples need to be merged for further processing to be
    identified as one single techincal rep.

        cd /media/disk6/ARF/RNASEQ/rawreads/LGEPfasta/lane11
        ls -d S* > list.lane
        while read file; do cat ./$file/*.fastq.gz > $file.fastq.gz; done < list.lane

        cd /media/disk6/ARF/RNASEQ/rawreads/LGEPfasta/lane12
        ls -d S* > list.lane
        while read file; do cat ./$file/*.fastq.gz > $file.fastq.gz; done < list.lane

        ....
        ....

2. Run fastqc on raw files
--------------------------

-   Change into each of the 12 directories; create list of .fastq.gz
    files; make fastqc output dir and run fastqc

        cd /media/disk6/ARF/RNASEQ/rawreads/EWfasta/lane10
        ls *.fastq.gz > lane10.list
        mkdir fastqc
        while read file; do /media/disk6/ARF/RNASEQ/software/FastQC/fastqc $file -o ./fastqc; done < lane10.list

        cd /media/disk6/ARF/RNASEQ/rawreads/EWfasta/lane12
        ls *.fastq.gz > lane12.list
        mkdir fastqc
        while read file; do /media/disk6/ARF/RNASEQ/software/FastQC/fastqc $file -o ./fastqc; done < lane12.list
        ...

-   Change into each of the 18 directories; create list of .fastq.gz
    files; make fastqc output dir and run fastqc

        cd /media/disk6/ARF/RNASEQ/rawreads/LGEPfasta/lane10
        ls *.fastq.gz > lane10.list
        mkdir fastqc
        while read file; do /media/disk6/ARF/RNASEQ/software/FastQC/fastqc $file -o ./fastqc; done < lane10.list


        cd /media/disk6/ARF/RNASEQ/rawreads/EWfasta/lane12
        ls *.fastq.gz > lane12.list
        mkdir fastqc
        while read file; do /media/disk6/ARF/RNASEQ/software/FastQC/fastqc $file -o ./fastqc; done < lane12.list
        ...

**Noticed there is rRNA contamination in Sample 71 in EW lane 09.**

3. Filter and trim raw reads
----------------------------

**rRNA contamination**

-   Accessed www.arb-silva.de and searched for Pinus taeda

-   805 Accessions showed up in results and were downloaded. The
    original .tgz is located at

        /media/disk6/ARF/RNASEQ/software/arb-silva.de_2018-01-28_id495190.tgz

-   The .fasta file inside the .tgz was placed in:

        /media/disk6/ARF/RNASEQ/software/bbmap/resources

-   bbDuk can be used to do overall quality trim, as well as trim first
    10 bases, remove adapters, and matches to rRNA. Test it out on EW
    sample 84 from lane10

        ./software/bbmap/bbduk.sh -Xmx2g in=./rawreads/EWfasta/lane10/84_S29_L002_R1_001.fastq.gz out=./8429l2.fq.gz      ref=./software/bbmap/resources/truseq.fa.gz,./software/bbmap/resources/rRNA_lob_pine.fasta forcetrimleft=10 t=1 ktrim=r k=25 mink=12 hdist=1 qtrim=rl trimq=20 minlength=50 maxlength=125

-   Initially Sample 84 had rRNA and adapter contamination. After
    running bbduk there was still one overresprestned sequence which
    when blasted hit to a tumor response protein. Since this is
    biologicaly relevenat it should be OK to stay. We can use bbduk now
    to process all raw reads.

-   Article suggest low to minimium quality trimming works best:
    <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4766705/>

-   Also minimum read length for DE expression roughly 50 works:
    <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4531809/>

-   First run on the EW technical reps a single sample at a time

        cd /media/disk6/ARF/RNASEQ
        while read file; do ./software/bbmap/bbduk.sh -Xmx2g in=./rawreads/EWfasta/lane02/$file     out=./trimmedfiltreads/EW/lane02/$file ref=./software/bbmap/resources/truseq.fa.gz,./software/bbmap/resources/rRNA_lob_pine.fasta forcetrimleft=10 t=5 ktrim=r k=25 mink=12 hdist=1 hdist2=0 qtrim=rl trimq=20 minlength=50 ; done < ./rawreads/EWfasta/lane02/lane02.list

        while read file; do ./software/bbmap/bbduk.sh -Xmx3g in=./rawreads/EWfasta/lane02/$file out=./trimmedfiltreads/EW/lane02/$file ref=./software/bbmap/resources/truseq.fa.gz,./software/bbmap/resources/rRNA_lob_pine.fasta forcetrimleft=10 t=5 ktrim=r k=25 mink=12 hdist=1 hdist2=0 qtrim=rl trimq=20 minlength=50 ; done < ./rawreads/EWfasta/lane02/lane02.list

-   Output is located at:

        /media/disk6/ARF/RNASEQ/trimmedfiltreads/EW/lane01 

-   Now run on the LGEP technical reps a single sample at a time

        while read file; do ./software/bbmap/bbduk.sh -Xmx3g in=./rawreads/LGEPfasta/lane01/$file out=./trimmedfiltreads/LGEP/lane01/$file ref=./software/bbmap/resources/truseq.fa.gz,./software/bbmap/resources/rRNA_lob_pine.fasta forcetrimleft=10 t=5 ktrim=r k=25 mink=12 hdist=1 hdist2=0 qtrim=rl trimq=20 minlength=50 ; done < ./rawreads/LGEPfasta/lane01/lane01.list


        while read file; do ./software/bbmap/bbduk.sh -Xmx3g in=./rawreads/LGEPfasta/lane09/$file out=./trimmedfiltreads/LGEP/lane09/$file ref=./software/bbmap/resources/truseq.fa.gz,./software/bbmap/resources/rRNA_lob_pine.fasta forcetrimleft=10 t=5 ktrim=r k=25 mink=12 hdist=1 hdist2=0 qtrim=rl trimq=20 minlength=50 ; done < ./rawreads/LGEPfasta/lane09/lane09.list

-   Multiple lanes ran one after the other:

        while read file; do ./software/bbmap/bbduk.sh -Xmx3g in=./rawreads/LGEPfasta/lane07/$file out=./trimmedfiltreads/LGEP/lane07/$file ref=./software/bbmap/resources/truseq.fa.gz,./software/bbmap/resources/rRNA_lob_pine.fasta forcetrimleft=10 t=5 ktrim=r k=25 mink=12 hdist=1 hdist2=0 qtrim=rl trimq=20 minlength=50 ; done < ./rawreads/LGEPfasta/lane07/lane07.list; while read file; do ./software/bbmap/bbduk.sh -Xmx3g in=./rawreads/LGEPfasta/lane08/$file out=./trimmedfiltreads/LGEP/lane08/$file ref=./software/bbmap/resources/truseq.fa.gz,./software/bbmap/resources/rRNA_lob_pine.fasta forcetrimleft=10 t=5 ktrim=r k=25 mink=12 hdist=1 hdist2=0 qtrim=rl trimq=20 minlength=50 ; done < ./rawreads/LGEPfasta/lane08/lane08.list

All lanes are processed the same way by changing the lane numbers in the
bbduk command.

Note: After processing Sample\_101b in lane10 for LGEP file was much
smaller, re-ran again on that file at a later time so file date is
newer. Same with Sample\_57h in lane13 for LGEP, except issue was that
fq.gz file was not complete.

4. FastQC can be re-done on filtered and trimmed files
------------------------------------------------------

5. Create complete experimental design matrix
---------------------------------------------

### Check sample indicies

Make link to: references/exptdesign/sequencing/

### Integrate phenotypes, pedigree, and sequencing data

Make link to: references/exptdesing/sequencing/

6. Align technical replicates with Salmon
-----------------------------------------

\*First index the 86k transcriptome:

\*p = 5 (five threads)

\*type = quasi (quasi-based mapping recommended)

\*k = 25 (A k of 31 is good for minimum read length = 75, but our
minimum read length is 50)

\*perfectHash (takes a little longer to build index but faster for quant
command)

    ./software/Salmon/bin/salmon index -t ./references/Pita.86ktxptome.fasta -i ./references/Pita.86k -p 5 --perfectHash --type "quasi" -k 25

\*p = 40 (use 40 threads)

\*l = U (library is unstranded)

\*numBootstraps = 10

    Salmon has the ability to optionally compute bootstrapped abundance estimates. This is done by resampling (with replacement) from the counts assigned to the fragment equivalence classes, and then re-running the optimization procedure, either the EM or VBEM, for each such sample. The values of these different bootstraps allows us to assess technical variance in the main abundance estimates we produce. Such estimates can be useful for downstream (e.g. differential expression) tools that can make use of such uncertainty estimates. This option takes a positive integer that dictates the number of bootstrap samples to compute. The more samples computed, the better the estimates of varaiance, but the more computation (and time) required.

\*seqBias

    Passing the --seqBias flag to Salmon will enable it to learn and correct for sequence-specific biases in the input data. Specifically, this model will attempt to correct for random hexamer priming bias, which results in the preferential sequencing of fragments starting with certain nucleotide motifs. By default, Salmon learns the sequence-specific bias parameters using 1,000,000 reads from the beginning of the input. If you wish to change the number of samples from which the model is learned, you can use the --numBiasSamples parameter. Salmon uses a variable-length Markov Model (VLMM) to model the sequence specific biases at both the 5’ and 3’ end of sequenced fragments. This methodology generally follows  that of Roberts et al. [2], though some details of the VLMM differ.

\*Execute Salmon in R

The experimental data set previously created can be used to run Salmon
on the techincal replicates.

    # Load experimental meta data and set working directory
    load("/media/disk6/ARF/RNASEQ/Breeding-Value-Prediction/disk6directory/resources/expt.dat.720.RData")
    setwd("/media/disk6/ARF/RNASEQ")

The first 288 rows correspond to the EW samples and the remaining 432
rows correspond to the LGEP samples. For each of the EW replicates, we
can use the experimental data to paste together a command which will
execute the alignment of a single technical replicate.

    for(each.file in 1:288){
      file <- paste0("/media/disk6/ARF/RNASEQ/trimmedfiltreads/86k/EW/",exp.dat$salmon_input_path[each.file])
      out.name <- paste0("Sample_",exp.dat$sample_id[each.file],"_",exp.dat$lane[each.file])
      true.lane <- gsub(pattern = "ew_",replacement = "",x = exp.dat$lane[each.file])
      cmd <- paste0("./software/Salmon/bin/salmon quant --seqBias -r ",file, " -i ./references/Pita.86k -o ./counts/86kSalmon/EW/", true.lane,"/",out.name," -p 40 -l U --numBootstraps 10 --fldMean ", exp.dat$fl_mean[each.file]," --fldSD 21")
      system(cmd)
    }  

The same can be done for the remaining rows which correspond to the
LGEP.

    for(each.file in 289:720){
      file <- paste0("/media/disk6/ARF/RNASEQ/trimmedfiltreads/86k/LGEP/",exp.dat$salmon_input_path[each.file])
      out.name <- paste0("Sample_",exp.dat$sample_id[each.file],"_",exp.dat$lane[each.file])
      cmd <- paste0("./software/Salmon/bin/salmon quant --seqBias -r ",file, " -i ./references/Pita.86k -o ./counts/86kSalmon/LGEP/", exp.dat$lane[each.file],"/",out.name," -p 40 -l U --numBootstraps 10 --fldMean ", exp.dat$fl_mean[each.file]," --fldSD 21")
      system(cmd)
    }  

-   All output files are located at:

        EW: /media/disk6/ARF/RNASEQ/counts/86kSalmon/EW/
        LGEP: /media/disk6/ARF/RNASEQ/counts/86kSalmon/LGEP/

7. Align biological replicates with Salmon
------------------------------------------

Similar to how the technical replicates were aligned, the experimental
meta data can be used to pass the 3 techincal replciates of each
biological replicate to Salmon and output as the biological replicate
counts. The mean fragment insert size was used across technical
replicates.

-   Load experimental meta-data and set working directory

        load("/media/disk6/ARF/RNASEQ/Breeding-Value-Prediction/disk6directory/resources/expt.dat.720.RData")
        setwd("/media/disk6/ARF/RNASEQ")

-   Identify the unique set of animal id names (biological replicates).

        ew.animals <- unique(as.character(expt.dat.720$animal_id[which(expt.dat.720$batch == "EW")]))

-   Now for each of the 96 EW biological replicates:

    -   identify the rows which correspond to the single biological
        replicate

    -   paste the input path of the exp-data file of the 3 technical
        replicates to the remaining full path

    -   paste the 3 full path techincal replicates together

    -   specify the file output name as Sample\_&lt;animal.name&gt;

    -   calculate the mean FL among the technical replicates

    -   plug information into system command and execute Salmon

<!-- -->

    for(each.file in 1:96){
     select.rows <-  which(as.character(expt.dat.720$animal_id) %in% ew.animals[each.file])

      file.cat.names <- (unlist( lapply(1:length(select.rows),function(x) { paste0("/media/disk6/ARF/RNASEQ/trimmedfiltreads/86k/EW/",expt.dat.720$salmon_input_path[select.rows[x]])}) ))
      
      file.cat.names <- paste(file.cat.names[1],file.cat.names[2],file.cat.names[3],sep = " ")
      
      out.name <- paste0("Sample_",ew.animals[each.file])
     FL.value <- mean(expt.dat.720$fl_mean[select.rows])
      cmd <- paste0("./software/Salmon/bin/salmon quant --seqBias -r ",file.cat.names, " -i ./references/Pita.86k -o ./counts/86kSalmon/bio_EW/",out.name," -p 50 -l U --numBootstraps 10 --fldMean ", FL.value," --fldSD 21")
      system(cmd)
    }  

-   The same is done for the LGEP

<!-- -->

    lgep.animals <- unique(as.character(expt.dat.720$animal_id[which(expt.dat.720$batch == "LGEP")]))
    for(each.file in 1:144){
     select.rows <-  which(as.character(expt.dat.720$animal_id) %in% lgep.animals[each.file])

      file.cat.names <- (unlist( lapply(1:length(select.rows),function(x) { paste0("/media/disk6/ARF/RNASEQ/trimmedfiltreads/86k/LGEP/",expt.dat.720$salmon_input_path[select.rows[x]])}) ))
      
      file.cat.names <- paste(file.cat.names[1],file.cat.names[2],file.cat.names[3],sep = " ")
      
      out.name <- paste0("Sample_",lgep.animals[each.file])
     FL.value <- mean(expt.dat.720$fl_mean[select.rows])
      cmd <- paste0("./software/Salmon/bin/salmon quant --seqBias -r ",file.cat.names, " -i ./references/Pita.86k -o ./counts/86kSalmon/bio_LGEP/",out.name," -p 55 -l U --numBootstraps 10 --fldMean ", FL.value," --fldSD 21")
      system(cmd)
    }  

-   Output files are located at:

<!-- -->

    EW: /media/disk6/ARF/RNASEQ/counts/86kSalmon/bio_EW/
    LGEP: /media/disk6/ARF/RNASEQ/counts/86kSalmon/bio_LGEP/
