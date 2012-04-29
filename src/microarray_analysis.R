# Change to your respective trunk
TRUNK.DIR<-"~/cap-6545-microarray-analysis/"

COLOR.PALETTE = c("green","red","blue","orange","purple","pink")

GetFilePaths <-function()
{
	# A bunch of static data for our datasets that we will be reading in
	my.expressionProfiles <- c("data/COLON/formatted/expression_profiles.csv",
	    					   "data/LEUKEMIA/formatted/expression_profiles.csv",
						       "data/LYMPHOMA/formatted/expression_profiles.csv",
							   "data/PROSTATE/formatted/expression_profiles.csv"
							   )
	my.classifications <- c("data/COLON/formatted/classification.txt",
							"data/LEUKEMIA/formatted/classification.txt",
							"data/LYMPHOMA/formatted/classification.txt",
							"data/PROSTATE/formatted/classification.txt"
							)
	my.genes <- c("data/COLON/formatted/genes.txt",
				  "data/LEUKEMIA/formatted/genes.txt",
				  "data/LYMPHOMA/formatted/genes.txt",
				  "data/PROSTATE/formatted/genes.txt"
				 )                    
				
	data.frame( expressionProfiles = my.expressionProfiles,
				classifications = my.classifications,
                genes = my.genes)
}   

GetData <-function()
{
	my.optimizedSigmas <- c(0.14,
					 	    0.10
						    0.001,
						    0.090
						    )	
					
	my.names <- c("Colon Dataset", 
				  "Leukemia Dataset",
				  "Lymphoma Dataset",
				  "Prostate Dataset"
				 )

	# NOTE: Not used yet
	my.kpcaAccuracies <- c(0.0,
						   0.0,
						   0.0,
						   0.0
						  )
	# NOTE: Not used yet
	my.biplotAccuracies <- c(0.0,
		                     0.0,
							 0.0,
							 0.0
							)						
						
	
	list( optimizedSigmas = my.optimizedSigmas,
		  names = my.names,
          kpcaAccuracies = my.kpcaAccuracies,
		  biplotAccuracies = my.biplotAccuracies)
}

				
# Sigh, R does not have good support for telling the path of the file that
# you're executing
setwd(TRUNK.DIR)

# We need the function GE.plot from the ge-biplot library.
source("public/ge-biplot.r", chdir=TRUE)
source("src/kpca_biplot.R", chdir=TRUE)

options.usepdf = TRUE


readAll<-function()
{
	my.expressionProfiles<-list()
	my.classifications<-list()
	my.colors<-list()
	for (i in 1:nrow(DATASET_FILES))
	{
		dataset.dirs <- DATASET_FILES[i,]
	    
		# Read in the expression profile
		expr.file <- as.matrix(dataset.dirs[which(names(dataset.dirs) == "expressionProfiles")])
		# Append matrix
		my.expressionProfiles[[i]] <- as.matrix(read.csv(file= expr.file,
														 header=FALSE,
														 sep=","))
	    
	    
		# Read in the labels for each expression profile
		class.file <- as.matrix(dataset.dirs[which(names(dataset.dirs) == "classifications")])
		# Append vector
		my.classification <- as.matrix(read.table(file= class.file,
								                  header=FALSE))[, 1]
		my.classifications[[i]] <- my.classification
		
		# Next, pair a color to each classification using a global color palette lookup table
		{
			class.factors <- factor(my.classification)
			class.factorslength <- length(unique(class.factors))
			my.color <- rep("", length(my.classification))
	        
			for (colorIndex in 1:class.factorslength)
			{
				my.color[my.classification == unique(class.factors)[colorIndex]]= COLOR.PALETTE[colorIndex]
			}
			my.colors[[i]] <- my.color
		}
	}

	list( expressionProfiles = my.expressionProfiles,
		  colors = my.colors,
		  classifications = my.classifications)
}


plotDifferentAlpha <- function(data, indexToOutputAlphas)
{
	for (currentAlpha in seq(from=0.0, to=1.0, by=0.10))
	{	
		# Perform our KPCA function
		kpca.data <- our.kpca(data$expressionProfiles[[indexToOutputAlphas]],
 							  kernel = "rbfdot",
							  kpar = list(sigma = data$optimizedSigmas[indexToOutputAlphas]), # Plot optimized value
 							  features = 2,
 							  alpha = currentAlpha)

				

		if (options.usepdf == TRUE)
		{
			kpca.pdfName <- paste(data$names[indexToOutputAlphas],"_alpha_",currentAlpha,"_kpca_biplot.pdf", sep="")
			pdf(kpca.pdfName)
		}
						
		# Plot KPCA GE-Biplot with optimized sigma
		kpca.biplot <- GE.plot(kpca.data$gene.expressions,
							   kpca.data$micro.arrays,
							   tit=   paste(data$names[indexToOutputAlphas],"KPCA-Biplot"),
							   clabs= data$classifications[[indexToOutputAlphas]],
							   glabs= "~",
							   cclr=  data$colors[[indexToOutputAlphas]],
							   gclr=  "black")      
		
		#DATASET.KPCA_ACCURACIES[i] <- doValidation(kpca.data$micro.arrays, data$classifications[[indexToOutputAlphas]])
		
		
		if (options.usepdf == TRUE)
		{
			dev.off()
		}
		
	}	
}



plotBest<-function(data)
{
	for (i in 1:length(data$names))
	{	
		# Perform our KPCA function
		kpca.data <- our.kpca(data$expressionProfiles[[i]],
 							  kernel = "rbfdot",
							  kpar = list(sigma = data$optimizedSigmas[i]), # Plot optimized value
 							  features = 2,
 							  alpha = 0.5)

				

		if (options.usepdf == TRUE)
		{
			kpca.pdfName <- paste(data$names[i],"_kpca_biplot.pdf")
			pdf(kpca.pdfName)
		}
						
		# Plot KPCA GE-Biplot with optimized sigma
		kpca.biplot <- GE.plot(kpca.data$gene.expressions,
							   kpca.data$micro.arrays,
							   tit=   "KPCA-Biplot",
							   clabs= data$classifications[[i]],
							   glabs= "~",
							   cclr=  data$colors[[i]],
							   gclr=  "black")      
		
		#DATASET.KPCA_ACCURACIES[i] <- doValidation(kpca.data$micro.arrays, data$classifications[[i]])
		
		
		if (options.usepdf == TRUE)
		{
			dev.off()
		}
						
		if (options.usepdf == TRUE)
		{
			biplot.pdfName <- paste(data$names[i],"_biplot.pdf")
			pdf(biplot.pdfName)
		}
						
		# Plot PCA GE-Biplot
		normal.biplot <- GE.biplot(data$expressionProfiles[[i]],
								   tit="GE-Biplot: Tumor Data",
								   clabs=data$classifications[[i]],
								   cclr=data$colors[[i]],
								   gclr="black",
								   opt.log=T,
								   opt.stand=T,
								   cex.lab=.8)
		
		#DATASET.BIPLOT_ACCURACIES[i] <- doValidation(normal.biplot$Chips, dataset.classifications)
			
		if (options.usepdf == TRUE)
		{
			dev.off()
		}
	}
}

doValidation <- function(microarrayMatrix,classifications,name,sigmaVal)
{	
	runningAccuracy <- 0
	numOfRuns       <- 2
	# Repeat 10 times to converge to the actual prediction
	for (runNumber in 1:numOfRuns)
	{	
		# 20-fold validation, auto-tuned parameters
		kvsmModel <- ksvm(microarrayMatrix, classifications, type="C-svc", C=1, cross=20)
		
		classificationPredictions <- predict(kvsmModel, microarrayMatrix)
		
		# How many did it correctly predict?
		currentAccuracy <- sum(classificationPredictions==classifications)/length(classifications)
		runningAccuracy <- runningAccuracy + currentAccuracy
	}
    library(stringr)
	png(file=paste("graphs/",name,"_",str_replace(sigmaVal,"\\.","_"),".png",sep=""))
    #plot the last one just to see what a random run would potentially look like
	plot(kvsmModel, data=microarrayMatrix)
    dev.off()
        
	accuracy <- runningAccuracy / numOfRuns
	
	accuracy
}

optimizeDataSet <- function(data,indexToOptimize)
{
	sigmaVector     <- c(0.001, 0.010, 0.050, 0.100, # low values
					 0.075, 0.080, 0.090, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, # normal values
					 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00) # high values

	optimizeSigma(data$expressionProfiles[[indexToOptimize]],
				  data$classifications[[indexToOptimize]],
				  sigmaVector,
                  data$names[indexToOptimize])
}

optimizeSigma <- function(expressionProfile,classifications,sigmaVector,name="default")
{
	accuracyVector  <- rep(0,length(sigmaVector))
	
	sigmaAccuracies <- list(sigmas = sigmaVector, accuracies = accuracyVector)
	
	for (i in 1:length(sigmaAccuracies$sigmas))
	{
		kpca.data <- our.kpca(expressionProfile,
						      kernel = "rbfdot",
 						      kpar = list(sigma = sigmaAccuracies$sigmas[i]),
 						      features = 2,
 							  alpha = 0.5)
		sigmaAccuracies$sigmas[i]
		sigmaAccuracies$accuracies[i] <- doValidation(kpca.data$micro.arrays, classifications, name,sigmaAccuracies$sigmas[i])
	}
	
	sigmaAccuracies
}


# Make a nice little dataframe out of all the paths to our datasets
DATASET_FILES <- GetFilePaths()
DATASET_DATA  <- GetData() 

# Always read all
my.data <- readAll()
DATASET_DATA$expressionProfiles <- my.data$expressionProfiles
DATASET_DATA$colors             <- my.data$colors
DATASET_DATA$classifications    <- my.data$classifications

#best.data <- plotBest(DATASET_DATA)

optimizeDataSet(DATASET_DATA,1)