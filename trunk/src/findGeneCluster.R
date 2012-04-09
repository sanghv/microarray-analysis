#
#	This finds the number of gene clusters for a given data set
#	This code is ugly. It was written quickly.
#
#	USAGE: In R Interactive Session:
#		source("findGeneCluster.R")
#	
#		# used to get predicted clusters
#		predictCluster(sigmvalue), see different values of sigma on a
#			a scatter plot in kernel space. This uses an information theoretic approach.
#
#		# used to get a biplot to see how it looks
#		makeBiPlot()
#
#		# used to find the best parameters for the kernel
#		findBestSplit()
#
#######################################################################

# information theoretic measure for clustering
source("jump.R")

# Change to your respective trunk
TRUNK.DIR<-"~/cap-6545-microarray-analysis/"
setwd(TRUNK.DIR)

# We need the function GE.plot from the ge-biplot library.
source("public/ge-biplot.r", chdir=TRUE)
source("src/kpca_biplot.R", chdir=TRUE)

COLOR.PALETTE = c("green","red","blue","orange","purple","pink")

# A bunch of static data for our datasets that we will be reading in
DATA.EXPRESSIONPROFILES <- c("data/COLON/formatted/expression_profiles.csv") 
DATA.CLASSIFICATION <- c("data/COLON/formatted/classification.txt")
DATA.GENES <- c("data/COLON/formatted/genes.txt")
                    
# Make a nice little dataframe out of all the paths to our datasets
DATASETS <- data.frame( expressionProfiles = DATA.EXPRESSIONPROFILES,
                        classifications = DATA.CLASSIFICATION,
                        genes = DATA.GENES)

#################################################################
# used for setting stuff up
#################################################################

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
			
	# Perform our KPCA function
	#kpca.data <- our.kpca(dataset.expressionProfiles, kernel = "rbfdot", kpar = list(sigma = 0.1), features = 2, alpha = 0.5)
	kpca.data <- 0
}

##################################################################
# Helper function ignore. Used for kmeans analysis
##################################################################
"km" <- function(x, output){
	y = jump(x)
	kms <- kmeans(x, y)
	print("Number of Clusters")
	print(y)
	
	if (output)
		plot(x, col=kms$cluster)	

}

###################################################################
# Make a nice graph of scatter points and clusters
###################################################################
"predictClusters" <- function(sigV, output=T) {
	
	kpca.data <- our.kpca(dataset.expressionProfiles, kernel = "rbfdot", kpar = list(sigma = sigV), features = 2, alpha = 0.5)
	y <- kpca.data$gene.expressions
	km(y, output)
	
	kpca.data
}

####################################################################
# make a ge-biplot to view what the sigma value gives
####################################################################
"makeBiPlot" <- function(kpca.data) {

	kpca.biplot <- GE.plot(kpca.data$gene.expressions,
                            kpca.data$micro.arrays,
                            tit=   "KPCA-Biplot",
                            clabs= dataset.classifications,
                            glabs= "~",
                            cclr=  dataset.classification.colors,
                            gclr=  "black")
}

#####################################################################
# SVM for finding best seperation 
#
######################################################################
"findBestSplit" <- function(i){
	
	y <- do.call(rbind, as.list(dataset.classifications))	
	
	#for (i in c(0.3)){
		dataSet <- predictClusters(i, output=F)
		x <- dataSet$gene.expressions	
		n <- nrow(x)	

		ntrain <- round(n*0.8)
		tindex <- sample(n, ntrain)

		xtrain <- x[tindex,]
		xtest <- x[-tindex,]	

		ytrain <- y[tindex]	
		ytest <- y[-tindex,]	
		
		supportVM <- ksvm(xtrain, ytrain, type="C-svc", kernel='vanilladot', C=100, scaled=c())
		
		ypred = predict(supportVM, xtest)

		print("Accuracy")
		print(sum(ypred==ytest)/length(ytest))

	#}
}
