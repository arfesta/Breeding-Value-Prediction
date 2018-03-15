# Step 3: Exploratory Data Analysis of normalized reads

  - Given that we will want to use the normalzied counts for the purposes of feature selection and ultimately prediction, we can start with taking a look at some summary stats and plots of how the normalized data sets return if they are consolidated to family mean or biological mean, as well as get aglimpse into the affect of applying any transformation proceducre.
  
  
  - To accomplish this, we will first alternate the use of two functions: One being [```descdist```](https://www.rdocumentation.org/packages/fitdistrplus/versions/1.0-8/topics/descdist) from the [```fitdistrplus```]
  9https://cran.r-project.org/web/packages/fitdistrplus/fitdistrplus.pdf) R package and the other being ```hist()``` from the base package.  The ```descdist``` function will provide a [skewness-kurtotois plot](http://www.itl.nist.gov/div898/handbook/eda/section3/eda35b.htm)
