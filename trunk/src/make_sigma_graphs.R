
# Get some decent anti-aliasing  up in this graph. 
library('Cairo')

library(ggplot2)

makeSigmaGraph<-function(plotTitle,plotDataset)
{
    #png(file=,sep=""))
    ggplot(data=plotDataset, aes(x=sigmas, y=acc, group=1)) +
                          geom_line(size=1.0) +
                          geom_smooth(method=loess, size=1.0, span=0.40) +
                          opts(title=paste("Optimized Sigma using KSVM for",plotTitle)) +
                          labs(x="Sigma Value", y="Prediction Accuracy")
    ggsave(filename=paste("../graphs/sigma_graph_for_",plotTitle,".png", sep=""))    
}

###################
# LYMPHOMA DATASET
#
# MOST APPROPRIATE SIGMA = 0.001
###################
makeSigmaGraph("Lymphoma",data.frame(sigmas = c( 0.001, 0.010, 0.050, 0.100, 0.075, 0.080, 0.090, 0.100, 0.110, 0.120, 0.130, 0.140,
								   0.150, 0.200, 0.300, 0.400, 0.500, 0.600, 0.700, 0.800, 0.900, 1.000),
					                 acc = c( 0.8831169, 0.8701299, 0.8506494, 0.8376623, 0.8376623, 0.8441558, 0.8441558,
								              0.8441558, 0.8376623, 0.8311688, 0.8246753, 0.8311688, 0.8311688, 0.8441558,
								              0.8766234, 0.8571429, 0.9090909, 0.9090909, 0.9090909, 0.8896104, 0.8961039,
								              0.9090909)))

###################
# PROSTATE DATASET
#
# MOST APPROPRIATE SIGMA = 0.090
###################                                   
makeSigmaGraph("Prostate",data.frame(sigmas = c( 0.001, 0.010, 0.050, 0.100, 0.075, 0.080, 0.090, 0.100, 0.110, 0.120, 0.130, 0.140,
								   0.150, 0.200, 0.300, 0.400, 0.500, 0.600, 0.700, 0.800, 0.900, 1.000),
					                 acc = c( 0.6764706, 0.7500000, 0.6960784, 0.6715686, 0.7401961, 0.6372549, 0.7500000,
                                              0.7107843, 0.6323529, 0.6715686, 0.6764706, 0.6274510, 0.6470588, 0.6421569,
                                              0.6568627, 0.7254902, 0.7500000, 0.7843137, 0.7892157, 0.8088235, 0.7990196,
                                              0.8382353)))
                                            

###################
# PROSTATE DATASET
#                                   
# MOST APPROPRIATE SIGMA = 0.14
##################
makeSigmaGraph("Tumor",data.frame(sigmas = c( 0.001, 0.010, 0.050, 0.100, 0.075, 0.080, 0.090, 0.100, 0.110, 0.120, 0.130, 0.140,
								   0.150, 0.200, 0.300, 0.400, 0.500, 0.600, 0.700, 0.800, 0.900, 1.000),
					                 acc = c( 0.8709677, 0.8709677, 0.8870968, 0.8709677, 0.8790323, 0.8709677, 0.8870968,
                                              0.8709677, 0.8709677, 0.8790323, 0.8870968, 0.9032258, 0.8870968, 0.7983871,
                                              0.7016129, 0.7177419, 0.7741935, 0.6612903, 0.7500000, 0.8145161, 0.7903226,
                                              0.7983871)))


###################
# LEUKEMIA DATASET
#
# MOST APPROPRIATE SIGMA = ???
###################                                          
#makeSigmaGraph("Leukemia",data.frame(sigmas = c(),
#					                 acc = c()))
        