# Change to your respective trunk
#TRUNK.DIR<-"C:/Users/King Charles/workspace/cap-6545-microarray-analysis/"
TRUNK.DIR<-"~/cap-6545-microarray-analysis/"

COLOR.PALETTE = c("green","red","blue","orange","purple","pink")

# A bunch of static data for our datasets that we will be reading in
DATA.EXPRESSIONPROFILES <- c("data/COLON/formatted/expression_profiles.csv", 
                             "data/LEUKEMIA/formatted/expression_profiles.csv",
                             "data/LYMPHOMA/formatted/expression_profiles.csv",
                             "data/PROSTATE/formatted/expression_profiles.csv")
DATA.CLASSIFICATION <- c("data/COLON/formatted/classification.txt",
                         "data/LEUKEMIA/formatted/classification.txt",
                         "data/LYMPHOMA/formatted/classification.txt",
                         "data/PROSTATE/formatted/classification.txt")
DATA.GENES <- c("data/COLON/formatted/genes.txt",
                "data/LEUKEMIA/formatted/genes.txt",
                "data/LYMPHOMA/formatted/genes.txt",
                "data/PROSTATE/formatted/genes.txt")                    
                    
# Make a nice little dataframe out of all the paths to our datasets
DATASETS <- data.frame( expressionProfiles = DATA.EXPRESSIONPROFILES,
                        classifications = DATA.CLASSIFICATION,
                        genes = DATA.GENES)


# Sigh, R does not have good support for telling the path of the file that
# you're executing
setwd(TRUNK.DIR)

# We need the function GE.plot from the ge-biplot library.
source("public/ge-biplot.r", chdir=TRUE)
source("src/kpca_biplot.R", chdir=TRUE)


for (i in 1:nrow(DATASETS))
{
    dataset.dirs <- DATASETS[i,]
    
    expr.file <- as.matrix(dataset.dirs[which(names(dataset.dirs) == "expressionProfiles")])
    dataset.expressionProfiles <- as.matrix(read.csv(file= expr.file,
                                                     header=FALSE,
                                                     sep=","))
    
    
    # Read in the labels for each expression profile
    class.file <- as.matrix(dataset.dirs[which(names(dataset.dirs) == "classifications")])
    dataset.classifications <- as.matrix(read.table(file= class.file,
                                                    header=FALSE))[, 1]
    # Next, pair a color to each classification using a global color palette lookup table
    {
        class.factors <- factor(dataset.classifications)
        class.factorslength <- length(unique(class.factors))
        dataset.classification.colors <- rep("", length(dataset.classifications))
        
        for (colorIndex in 1:class.factorslength)
        {
            dataset.classification.colors[dataset.classifications == unique(class.factors)[colorIndex]]= COLOR.PALETTE[colorIndex]
        }
    }
        
                               
    # TODO: Read in gene metadata. It's parsed and waiting for us.
    # ...
    # ...

    # Perform our KPCA function
    kpca.data <- our.kpca(dataset.expressionProfiles, kernel = "rbfdot",
                          kpar = list(sigma = 0.1), features = 2, alpha = 0.5)

    # Plot it 
    kpca.biplot <- GE.plot(kpca.data$gene.expressions,
                           kpca.data$micro.arrays,
                           tit=   "KPCA-Biplot",
                           clabs= dataset.classifications,
                           glabs= "~",
                           cclr=  dataset.classification.colors,
                           gclr=  "black")           
					
	# TODO: Add in regular GE Biplot, clean up code, etc. 
    
}

                 
                 
######################################################################
# OLD CODE BELOW, DIDN'T HAVE TIME TO EXTRACT WHAT I WANTED FROM IT  #
######################################################################
## Relative to your trunk directory
#I2000.PATH   <- "data/colon tumor/raw formatted data/I2000.csv"
#TISSUES.PATH <- "data/colon tumor/raw formatted data/tissues.txt"
#
#
#
## Read in formatted tumor data
#tumor.data <- as.matrix(read.csv(file=I2000.PATH,header=FALSE,sep=","))
#
## Read in tissue data as a vector
#tissue.data <- as.matrix(read.table(file=TISSUES.PATH,header=FALSE))[, 1]
#
## All Green by default, then change '-' to Red
#tissue.col<-rep("green",length(tissue.data))
#tissue.col[tissue.data=="-"]="red"
#
#kpca.tumor.data <- our.kpca(tumor.data, kernel = "rbfdot",
                            #kpar = list(sigma = 0.1), features = 2, alpha = 0.5)
#
## kpca.tumor.data is a list of the gene expression data and the microarrays
## projected onto the principal components. 
#
#options.usepdf = FALSE
#
#
#if (options.usepdf == FALSE)
#{
	#par(mfrow=c(1,2),mar=c(2,2,2,2),ask=F)
#}
#
## The kpca biplot from Reverter, Vegas, and Sanchez.
#if (options.usepdf == TRUE)
#{
	#pdf("kpca-biplot.pdf")
#}
#kpca.biplot <- GE.plot(kpca.tumor.data$gene.expressions,
                       #kpca.tumor.data$micro.arrays,
                       #tit = "KPCA-Biplot: Tumor Data", clabs = tissue.data,
                       #glabs = "~", cclr =tissue.col, gclr="black")
#if (options.usepdf == TRUE)
#{
	#dev.off()
#}			
## The linear biplot from Pittelkow and Wilson
#if (options.usepdf == TRUE)
#{
	#pdf("normal-biplot.pdf")
#}
#normal.biplot <- GE.biplot(tumor.data,tit="GE-Biplot: Tumor Data",
                         #clabs=tissue.data,cclr=tissue.col,gclr="black",
                         #opt.log=T,opt.stand=T,cex.lab=.8)
#if (options.usepdf == TRUE)
#{
	#dev.off()
#}
