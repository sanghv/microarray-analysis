Preprocessing done according to the preprocessing section in "SuppInfo_CCv3.pdf" in doc/research papers/PROSTATE DataSet Papers.

It states:
"All expression files in a given experiment were scaled to a reference file (generally the
file found to have the median value of expression) based upon the mean average
difference for all genes present on the microarrays. The scaled files used in each
experiment will be available at www-genome.wi.mit.edu/MPR/Prostate.  All genes with
average differences below the minimum threshold of 10 were set at the minimum
threshold.  The maximum threshold was set at 16,000.  After thresholding, the relative
variation of expression for each gene was determined by dividing the maximum
expression for the gene among all samples (Max) by the minimum expression (Min)
(Max/Min).  The absolute variation in expression was determined by subtracting the Min
from the Max (Max-Min).  Filtering parameters of 5-fold change (Max/Min) and absolute
difference of 50 (Max-Min) were used for all subsequent analysis."