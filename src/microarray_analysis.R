# Change to your respective trunk
TRUNK.DIR<-"C:/Users/King Charles/workspace/cap-6545-microarray-analysis/"
#TRUNK.DIR<-"~/cap-6545-microarray-analysis/"

COLOR.PALETTE = c("green","red","blue","orange","purple","pink")

GetFilePaths <-function()
{
	# A bunch of static data for our datasets that we will be reading in
	my.expressionProfiles <- c("data/COLON/formatted/expression_profiles.csv"#, <--- remove this too lol
	    					   #"data/LEUKEMIA/formatted/expression_profiles.csv",
						       #"data/LYMPHOMA/formatted/expression_profiles.csv",
							   #"data/PROSTATE/formatted/expression_profiles.csv"
							   )
	my.classifications <- c("data/COLON/formatted/classification.txt"#, <--- remove this too lol
							#"data/LEUKEMIA/formatted/classification.txt",
							#"data/LYMPHOMA/formatted/classification.txt",
							#"data/PROSTATE/formatted/classification.txt"
							)
	my.genes <- c("data/COLON/formatted/genes.txt"#, <--- remove this too lol
				  #"data/LEUKEMIA/formatted/genes.txt",
				  #"data/LYMPHOMA/formatted/genes.txt",
				  #"data/PROSTATE/formatted/genes.txt"
				 )                    
				
	data.frame( expressionProfiles = my.expressionProfiles,
				classifications = my.classifications,
                genes = my.genes)
}   

GetData <-function()
{
	my.optimizedSigmas <- c(0.14#, <--- remove this too lol
					 	    #0.10
						    #0.001,
						    #0.090
						    )	
					
	my.names <- c("Colon Dataset"#, <--- remove this too lol
				  #"Leukemia Dataset",
				  #"Lymphoma Dataset",
				  #"Prostate Dataset"
				 )

	# NOTE: Not used yet
	my.kpcaAccuracies <- c(0.0#,
						   #0.0,
						   #0.0,
						   #0.0
						  )
	# NOTE: Not used yet
	my.biplotAccuracies <- c(0.0#,
		                     #0.0,
							 #0.0,
							 #0.0
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

# NOTE: Probably want to do the tumor one since we get decent looking output for it

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

doValidation <- function(microarrayMatrix,classifications)
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
	#plot the last one just to see what a random run would potentially look like
	#plot(kvsmModel, data=microarrayMatrix)

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
				  sigmaVector)
}

optimizeSigma <- function(expressionProfile,classifications,sigmaVector)
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
		sigmaAccuracies$accuracies[i] <- doValidation(kpca.data$micro.arrays, classifications)
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

# Make sure the index corresponds to whatever dataset you want to plot alphas for.
# I'm gonna get what different alphas look like for the tumor dataset
plotDifferentAlpha(DATASET_DATA, 1)

# LYMPHOMA (BEST SIGMA = 0.001)
# $sigmas
# [1] 0.001 0.010 0.050 0.100 0.075 0.080 0.090 0.100 0.110 0.120 0.130 0.140
#[13] 0.150 0.200 0.300 0.400 0.500 0.600 0.700 0.800 0.900 1.000
#
#$accuracies
# [1] 0.8831169 0.8701299 0.8506494 0.8376623 0.8376623 0.8441558 0.8441558
# [8] 0.8441558 0.8376623 0.8311688 0.8246753 0.8311688 0.8311688 0.8441558 <-- saturation, bad values
#[15] 0.8766234 0.8571429 0.9090909 0.9090909 0.9090909 0.8896104 0.8961039 <-- saturation, bad values
#[22] 0.9090909                                                             <-- saturation, bad values

# TUMOR (BEST SIGMA = 0.14)
# $sigmas
# [1] 0.001 0.010 0.050 0.100 0.075 0.080 0.090 0.100 0.110 0.120 0.130 0.140
#[13] 0.150 0.200 0.300 0.400 0.500 0.600 0.700 0.800 0.900 1.000 5.000
#
#$accuracies
# [1] 0.8709677 0.8709677 0.8870968 0.8709677 0.8790323 0.8709677 0.8870968
# [8] 0.8709677 0.8709677 0.8790323 0.8870968 0.9032258 0.8870968 0.7983871
#[15] 0.7016129 0.7177419 0.7741935 0.6612903 0.7500000 0.8145161 0.7903226
#[22] 0.7983871 0.8064516

# PROSTATE (BEST SIGMA = 0.090)
#$sigmas
# [1] 0.001 0.010 0.050 0.100 0.075 0.080 0.090 0.100 0.110 0.120 0.130 0.140
#[13] 0.150 0.200 0.300 0.400 0.500 0.600 0.700 0.800 0.900 1.000
#
#$accuracies
# [1] 0.6764706 0.7500000 0.6960784 0.6715686 0.7401961 0.6372549 0.7500000
# [8] 0.7107843 0.6323529 0.6715686 0.6764706 0.6274510 0.6470588 0.6421569
#[15] 0.6568627 0.7254902 0.7500000 0.7843137 0.7892157 0.8088235 0.7990196	<-- saturation, bad values
#[22] 0.8382353																<-- saturation, bad values

# LEUKEMEIA (BEST SIGMA = 0.100)
#$sigmas
# [1] 0.001 0.010 0.050 0.100 0.075 0.080 0.090 0.100 0.110 0.120 0.130 0.140
#[13] 0.150 0.200 0.300 0.400 0.500 0.600 0.700 0.800 0.900 1.000
#
#$accuracies
# [1] 0.7500000 0.7777778 0.7569444 0.8194444 0.7638889 0.8125000 0.7708333
# [8] 0.8263889 0.7847222 0.7916667 0.7777778 0.8125000 0.7847222 0.8055556
#[15] 0.8680556 0.8472222 0.8750000 0.8333333 0.8333333 0.8333333 0.8194444	<-- saturation, bad values
#[22] 0.8125000																<-- saturation, bad values

#optimizeDataSet(DATASET_DATA,1)