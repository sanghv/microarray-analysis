Preprocessing done according to the preprocessing section in "Shipp_et_al_Supplementary_Information_v5.pdf" in doc/research papers/LYMPHOMA DataSet Papers.

It states:
"Preprocessing and Re-scaling
The raw expression data as obtained from Affymetrix's GeneChip is re-scaled to
account for different chip intensities. Each column (sample) in the data set was
multiplied by 1/slope of a least squares linear fit of the sample vs. the reference (the
first sample in the data set). This linear fit is done using only genes that have
'Present' (P) calls in both the sample being re-scaled and the reference.  (The P calls
are calculated by Affymetrix’s GENECHIP software and each P call represents a
gene with RNA “Present” as determined by the average difference analysis of
expression measurements from a gene’s set of probes on the microarray.) The
sample chosen as reference is a typical one (i.e. one with the number of "P" calls
closer to the average over all samples in the data set).
A ceiling of 16,000 units was chosen for all experiments because it is at this level
that we observe fluorescence saturation of the scanner; values above this cannot be
reliably measured.  We set a lower threshold for the expression levels to 20 units to
minimize noise effects while avoiding missing any potentially informative marker
genes.
These numbers are Affymetrix’s scanner “average difference” units. After this
preprocessing, gene expression values were subjected to a variation filter that
excluded genes showing minimal variation across the samples being analyzed.  The
variation filter tests for a fold-change and absolute variation over samples
(comparing max/min and max-min with predefined values and excluding genes not
obeying both conditions).  For maximum/minimum fold variation, we excluded genes
with less than 3-fold variation and, for maximum-minimum absolute variation, we
excluded genes with less than 100 units absolute variation."

We should still do normalization and scaling though since these are not done to the dataset to my knowledge.