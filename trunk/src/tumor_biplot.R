# Noobs to R read:
# http://www.johndcook.com/R_language_for_programmers.html

# Change to your respective trunk
TRUNK_DIR<-"C:/Users/King Charles/workspace/cap-6545-microarray-analysis/"
# Relative to your trunk directory
I2000_PATH   <- "data/colon tumor/raw formatted data/I2000.csv"
TISSUES_PATH <- "data/colon tumor/raw formatted data/tissues.txt"

# Sigh, R does not have good support for telling the path of the file that you're executing
setwd(TRUNK_DIR)
# Use the original GE-Biplot library until you've got a good handle on it.
source("public/ge-biplot.r", chdir=TRUE)

# Read in formatted tumor data
tumor_data <- as.matrix(read.csv(file=I2000_PATH,header=FALSE,sep=","))

# We need to do a matrix transpose to get it into the format it wants:
# rows = genes, cols = samples
# tumor_data <- t(tumor_data)

# Read in tissue data as a vector
tissue_data <- as.matrix(read.table(file=TISSUES_PATH,header=FALSE))[, 1]

# All Green by default, then change '-' to Red
tissue_col<-rep("green",length(tissue_data))
tissue_col[tissue_data=="-"]="red"

# i) log2 transformation
tumor_data<-log2(tumor_data)
# ii) chip (row) standardising these values
# iii) Column (gene) centering at 0.
# http://en.wikipedia.org/wiki/Centering_matrix
tumor_data<-t(scale(t(tumor_data)))

# I'm not super sure what the syntax '~.' is but how I read this whole statement:
# for each row of the data.frame (data.frame being like a table data structure)
# take each row, then apply 'laplacian dot kernel' with sigma being equal to 0.35
# then take the higest valued eigenvalue features from the results.
kpca_tumor_data <- kpca(~., data = data.frame(tumor_data), kernel="laplacedot", kpar = list(sigma=0.35), features=2)

# kpca_tumor_data is a S4 class.
# it has a bunch of information about it and you have to access it. You access
# the matrix by rotated(kpca_tumor_data)

plot(rotated(kpca_tumor_data), col=tissue_col, xlab="1st Principal Component",ylab="2nd Principal Component")


normal_biplot<-GE.biplot(tumor_data,tit="GE-Biplot: Tumor Data",clabs=tissue_data,cclr=tissue_col,gclr="black",
opt.log=F,opt.stand=F,cex.lab=.8)
