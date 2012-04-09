########################################################################################
#
# Hacked up code implementation for Sugar and James Entropy approximation for number of
# clusters to use in K-Means methods. The theory dictates a y vector should be formed
# where a matrix nxp (n samples....p genes) will give y = p/2. This is an awesome algo.
# 
# Don't bother reading this. It will make your brain melt. No attempts were made to be 
# clear. NO TIME!!!!!!!!!!!!!!!!!!
#
########################################################################################

# jump algorithm
"jump" <- function(data=NULL){
	
	K=10
	y = ncol(data)/2
	rand=10
	fits=NULL
	B=0
	dist=NULL

	if (!is.null(data)){
    	# Compute the kmeans fit to the data
  		if (is.null(fits))
  			fits <- kmeans.rndstart(data,K,rand)
  
		if (is.null(y))
    		y <- dim(data)[2]/2
  
		K <- length(fits)
  		n <- length(fits[[2]]$cluster)
  		p <- dim(fits[[2]]$cent)[2]
  		dist <- rep(0,K)
  
		# Compute the distortion associated with the kmeans fit
  		for (k in 1:K)
      		dist[k] <- sum(fits[[k]]$within)/(n*p)}
  
		# Call the compute.jump 
  		jump.results <- compute.jump(dist,y)
  		jump.results$fits <- fits
  
		# Implement bootstrap routine
  		if (B>0 & !is.null(data)){
  			n <- nrow(data)
  			boot.results <- matrix(0,length(y),K)
  			bootdist <- rep(0,K)
  
			for (b in 1:B){
    			# Make bootstrap data
    			bootdata <- data[sample(1:n,replace=T),]
    			
				# Get kmeans fit to the bootstrap data
    			bootfits <- kmeans.rndstart(bootdata,K,rand)
    
				# Compute bootstrap distortion and maximum jumps
    			for (k in 1:K)
      				bootdist[k] <- sum(bootfits[[k]]$within)/(n*p)
    
				bootmaxjumps <- compute.jump(bootdist,y,printresults=F)$maxjump
    
				for (j in 1:length(y))
      				boot.results[j,bootmaxjumps[j]] <-  boot.results[j,bootmaxjumps[j]]+1
  			}
  
		# Calculate proportions of each number of clusters chosen
  		jump.results$boot.result <- round(boot.results/B,3)
  
		for (j in 1:length(y))
    		print(paste(jump.results$boot.result[j,jump.results$maxjump[j]]*100,"% of bootstrap iterations corresponding to ",jump.results$maxjump[j], "clusters with Y=",y[j]))
     	}
  
		jump.results[1][1]$maxjump
}

#################################################################
# Helper function ---> ignore this
#################################################################
"compute.jump" <- function(dist,y,printresults=T){
	K <- length(dist)
	numb.y <- length(y)
	numbclust <- rep(0,numb.y)
	transdist <- matrix(0,numb.y,K+1)
	jumps <- matrix(0,numb.y,K)
    
	for (i in 1:numb.y){
    	# Compute the transformed distortion
      	transdist[i,] <- c(0,dist^(-y[i]))
      
		# Compute the jumps in transformed distortion
      	jumps[i,] <- diff(transdist[i,])
      	
		# Compute the maximum jump
      	numbclust[i] <- order(-jumps[i,])[1]
      
      	# Debug Print
      	#if (printresults)
    	#	print(paste("The maximum jump occurred at ",numbclust[i], "clusters with Y=",y[i]))
    }
	list(maxjump=numbclust,dist=dist,transdist=transdist[,-1],jumps=jumps)}

#################################################################
# Helper function ---> ignore this
#################################################################
"kmeans.rndstart" <- function(x, K, rand = 10){
  
	n <- nrow(x)
	fits <- list(list(withinss=sum((t(x) - apply(x, 2, mean))^2)))
	iter.max <- 10
  
	# Run kmeans for 2 to K clusters
	for (k in 2:K){
		bestss <- 10^20
    
	# Run kmeans rand times with different centers each time
    for(i in 1:rand) {
    	m <- nrow(x)
    	marker <- T
    
	while(marker) {
    	centers <- x[sample(m, k),  ]
    
		if (m < k) 
        	stop("more cluster centers than data points")
    
		# Fortran function to calculate kmeans solution
    	Z <- .Fortran("kmns", as.double(x), as.integer(m), as.integer(ncol(x)), 
        	centers = as.double(centers), as.integer(k), c1 = integer(m), 
        	integer(m), nc = integer(k), double(k), double(k), integer(k), 
        	double(m), integer(k), integer(k), as.integer(iter.max), 
        	wss = double(k), ifault = as.integer(0), PACKAGE = "stats")
    
		if (Z$ifault==2)
      		warning("did not converge in iter.max iterations",  call. = FALSE)
    
		if (Z$ifault==0)
      	marker=F
	}
    
	centers <- matrix(Z$centers, k)
    dimnames(centers) <- list(1:k, dimnames(x)[[2]])
    
	# Keep solution corresponding to lowest distortion
    if(sum(Z$wss) < bestss) {
		finalkmeans <- list(cluster = Z$c1, centers = centers, withinss = Z$wss, 
		size = Z$nc)
		bestss <- sum(finalkmeans$within)
	}
	}
	fits <- c(fits,list(finalkmeans))}
	fits
}
