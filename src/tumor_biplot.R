# Change to your respective trunk
TRUNK.DIR<-"~/cap-6545-microarray-analysis/"
# Relative to your trunk directory
I2000.PATH   <- "data/colon tumor/raw formatted data/I2000.csv"
TISSUES.PATH <- "data/colon tumor/raw formatted data/tissues.txt"

# Sigh, R does not have good support for telling the path of the file that
# you're executing
setwd(TRUNK.DIR)

# We need the function GE.plot from the ge-biplot library.
source("public/ge-biplot.r", chdir=TRUE)
source("src/kpca_biplot.R")

# Read in formatted tumor data
tumor.data <- as.matrix(read.csv(file=I2000.PATH,header=FALSE,sep=","))

# Read in tissue data as a vector
tissue.data <- as.matrix(read.table(file=TISSUES.PATH,header=FALSE))[, 1]

# All Green by default, then change '-' to Red
tissue.col<-rep("green",length(tissue.data))
tissue.col[tissue.data=="-"]="red"

kpca.tumor.data <- our.kpca(tumor.data, kernel = "rbfdot",
                            kpar = list(sigma = 0.1), features = 2, alpha = 0.5)

# kpca.tumor.data is a list of the gene expression data and the microarrays
# projected onto the principal components. 

# The kpca biplot from Reverter, Vegas, and Sanchez.
# Doesn't work as expected yet.
pdf("kpca-biplot.pdf")
kpca.biplot <- GE.plot(kpca.tumor.data$gene.expressions,
                       kpca.tumor.data$micro.arrays,
                       tit = "KPCA-Biplot: Tumor Data", clabs = tissue.data,
                       glabs = "~", cclr =tissue.col, gclr="black")
dev.off()

# The linear biplot from Pittelkow and Wilson
pdf("normal-biplot.pdf")
normal.biplot <- GE.biplot(tumor.data,tit="GE-Biplot: Tumor Data",
                         clabs=tissue.data,cclr=tissue.col,gclr="black",
                         opt.log=T,opt.stand=T,cex.lab=.8)
dev.off()
