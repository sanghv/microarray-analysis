library(MASS)
#________________________________________________________ 
#
#  Function to produce the GE-biplot as described
#  in Pittelkow & Wilson (2003). 
#________________________________________________________ 
GE.biplot<-
function(X,tit="",glabs="~", clabs=rownames(X), gclr="red",cclr="blue", axes=c(1,2), opt.gf=T,opt.plot=T, opt.log=T, opt.stand=T,...)
{
# X is a  matrix, organized as chip by genes
   X<-as.matrix(X)
   k<-max(axes)
   if(k<2 | k > min(dim(X))) 
     {cat("GE:ERROR: axes parameter incorrect");return}
   if(opt.log) {if ( min(X) < 0 ) 
     {cat("GE:ERROR: Can't take logs of negative or zero values");return}
     X<-log2(X)}
   if(opt.stand) {X<-t(scale(t(X)))}     
    
   gf<-c(I1=NA,I2=NA)
   mn<-colMeans(X)
   X <- sweep(X, 2, mn)
   sc<-dim(X)[1]-1
   W<-rep(sqrt(sc),k)
   
   GEout<-GE.biplt(X,k,0,W)
   
   if(opt.gf)
      {gf<-GE.gf(GEout$svdx$d,k)}
#
    if(opt.plot)
      {if(opt.gf) 
      {tit<- paste(tit, "\nI1=",round(gf[[1]],2), "I2=",round(gf[[2]],2))}
       GE.plot(GEout$Genes[,axes],GEout$Chips[,axes],tit,glabs,clabs,gclr,cclr, xlab=paste("Axis",axes[1]),ylab=paste("Axis",axes[2]),...)}
    
list(Genes=GEout$Genes,Chips=GEout$Chips,gf=gf,sv=GEout$svdx$d)
}


#________________________________________________________ 
#
#  Function to produce the Chip-plot as described
#  in Pittelkow & Wilson (2003). 
#________________________________________________________ 

GE.ChipPlot<-
function(X, tit="",glabs="~", clabs=rownames(X), gclr="red",cclr="blue", axes=c(1,2), opt.gf=T,opt.plot=T, opt.log=T, opt.stand=T,...)
{
# X is a  matrix, organized as chip by genes
   X<-as.matrix(X)
   k<-max(axes)
   if(k<2 | k > min(dim(X))) 
     {cat("GE:ERROR: axes parameter incorrect");return}
   if(opt.log) {if ( min(X) < 0 ) 
     {cat("GE:ERROR: Can't take logs of negative or zero values");return}
     X<-log2(X)}
   if(opt.stand) {X<-t(scale(t(X)))}     
    
   gf<-c(I1=NA,I2=NA)
   mn<-colMeans(X)
   X <- sweep(X, 2, mn)
   W<-rep(1,k)
   
   GEout<-GE.biplt(X,k,1,W)
   
   if(opt.gf)
      {gf<-GE.gf(GEout$svdx$d,k)}
#
    if(opt.plot)
      {if(opt.gf) 
      {tit<- paste(tit, "\nI1=",round(gf[[1]],2), "I2=",round(gf[[2]],2))}
       GE.plot(GEout$Genes[,axes],GEout$Chips[,axes],tit,glabs,clabs,gclr,cclr, xlab=paste("Axis",axes[1]),ylab=paste("Axis",axes[2]),...)}
    
list(Genes=GEout$Genes,Chips=GEout$Chips,gf=gf,sv=GEout$svdx$d)
}

#________________________________________________________ 
#
#  Function to produce the Gene-plot, by calling a user defined
#  function  which plots an icon at the gene point.
#  Examples can be found in Pittelkow & Wilson (2003). 
#________________________________________________________  
 
GE.GenePlot<-
function(X, tit="",icon="GE.Icon", genedata,gclr="black",axes=c(1,2), opt.gf=T, opt.log=T, opt.stand=T,...)
{
# X is a  numerical matrix, organized as chip by genes
# geneinfo is a numerical axis of means to be used for icons
   X<-as.matrix(X)
   
   k<-max(axes)
   if(k<2 | k > min(dim(X))) 
     {cat("GE:ERROR: axes parameter incorrect");return}
   if(opt.log) {if ( min(X) < 0 ) 
     {cat("GE:ERROR: Can't take logs of negative or zero values");return}
     X<-log2(X)}
   if(opt.stand) {X<-t(scale(t(X)))}     
   
   
   gf<-c(I1=NA,I2=NA)
   mn<-colMeans(X)
   X <- sweep(X, 2, mn)
   sc<-dim(X)[1]-1
   W<-rep(sqrt(sc),k)
   
   GEout<-GE.biplt(X,k,0,W)
   
   
   if(opt.gf)
      {gf<-GE.gf(GEout$svdx$d,k)}
#
      tit<- paste(tit, "\nI1=",round(gf[[1]],2), "I2=",round(gf[[2]],2)) 
      do.call(icon,list(Genes=GEout$Genes[,c(axes[1],axes[2])],tit=tit,genedata=genedata,
      gclr=gclr,iconsize=.1,xlab=paste("Axis ",axes[1]), ylab=paste("Axis",axes[2])))
   
list(Genes=GEout$Genes,Chips=GEout$Chips,gf=gf,sv=GEout$svdx$d)
}

#________________________________________________________
#  Function to produce the coordinates from a SVD 
#
#  Returns the coordinates and the SVD results.
#  Notes: 
#  If large problem, replace matrix multiplication with sweep.
#  Usually this function would not be called directly.
#  The coordinates will only be meaningful if appropriate 
#  transformation on the gene expressions have been carried
#  out on X before calling GE.biplt. 
#________________________________________________________ 
#
GE.biplt<-function(X,k=2,alpha=0,W=rep(sqrt(dim(X)[1]),k))
{
#  k dimensional solution
#  X is a  matrix, organized as chip by genes
   svdx<-La.svd(X,nu=k,nv=k)
   if (alpha==0)
      {C<-svdx$u%*%diag(W)
       G<-t(svdx$vt)%*%diag(svdx$d[1:k]/W)}
   if (alpha==1)
      {C<-svdx$u%*%diag(svdx$d[1:k]*W)
       G<-t(svdx$vt)%*%(diag(1/W))}
   list(Genes=G,Chips=C,svdx=svdx)
}
#_____________________________________________________ 

# Function to compute the Goodness of Fit statistics from 
# singular values returned by a call to GE.biplot for example.
#________________________________________________________ 

GE.gf<-function(sv,k)
{
d<-sv
i1<-(d^2)/sum(d^2)
i2<-(d^4)/sum(d^4)
c(I1=sum(i1[1:k]),I2=sum(i2[1:k]))
}
#________________________________________________________ 
#
# Function to plot the output from GE.biplt or GE.biplot
# Plots either text or R plot symbols 
# 
#________________________________________________________ 
GE.plot<-function(Genes,Chips,tit,glabs,clabs,gclr,cclr,xlab="Axis 1", ylab="Axis 2",...)
{
eqscplot(c(Genes[,1],Chips[,1]),c(Genes[,2],Chips[,2]),tol=.2,main=tit,type="n", xlab=xlab,ylab=ylab,...)
if(is.character(glabs[1]))
   text(Genes, labels=glabs, col=gclr,...)
if(is.numeric(glabs[1]))
   points(Genes, pch=glabs, col=gclr,...)
if(is.character(clabs[1]))
text(Chips[,1],Chips[,2],labels=clabs, col=cclr,...)
if(is.numeric(clabs[1]))
points(Chips[,1],Chips[,2],pch=clabs, col=cclr,...)
segments(0,0,Chips[,1],Chips[,2],col=cclr,lty="dotted",...)
abline(h=0,lty="dotted")
abline(v=0,lty="dotted")
}
#
#________________________________________________________ 
#
# Function to plot  comparative line segments for summary 
#     statistics
#________________________________________________________ 
GE.Icon<-function(Genes,tit,genedata,gclr="black",iconsize=.1,xlab="Axis 1", ylab="Axis 2",...)
{
  eqscplot(Genes,tol=.2,main=tit,type="n", xlab=xlab,ylab=ylab,...)
  ngrps<-dim(genedata)[2]
  iconbox<-(max(Genes)-min(Genes))*iconsize
  iconsc<-(max(genedata)-min(genedata))
  coords<-Genes
  for (gps in c(1:(ngrps-1))){
  rise<-(genedata[,gps+1]-genedata[,gps])*iconbox/iconsc
 coords2<-cbind((coords[,1]+iconbox/(ngrps-1)),(coords[,2]+rise)) 
segments(coords[,1],coords[,2],coords2[,1],coords2[,2],col=gclr,...)
coords<-coords2                           }
}
#
#________________________________________________________ 
#
# Function to identify genes on the current biplot
# Prints the label on the diagram  
# Returns the row position of the coordinates and the 
# columns of the gene expression matrix in the call to the 
# biplot function.
#________________________________________________________ 
GE.identifyGenes<-function(Genes,glabs=rownames(Genes),...)
{identify(Genes,labels=glabs,...)}  
#________________________________________________________ 
#
# Function to demonstrate the use of GE package  on the
# Golub data
#________________________________________________________ 
GE.demo<-function()
{
library(MASS) #Required for GE functions
library(multtest) 
data(golub)  
# NOTES on the data set
#   The Data is from the library multtest.
#   Only the classification into two tumor types is recorded.
#   To find more info on this set type 'help(golub)'.
#   The data appears to be  already chip centered and logged. 
#   I understand that a thresholding operation and a 
#   gene filtering step has been performed prior to 
#   chip standardizing and logging.
#   Results are different from Pittelkow and Wilson (1999) 
#   because of slightly different preprocessing.
#   The demo assumes the order of samples in the data set is 
#   the same as the WEB data.
#   Set up labels and colours for the chips
cells<-c("B-cell","T-cell","T-cell","B-cell","B-cell","T-cell",
"B-cell","B-cell","T-cell","T-cell","T-cell","B-cell","B-cell",
"T-cell","B-cell","B-cell","B-cell","B-cell","B-cell","B-cell",
"B-cell","B-cell","T-cell","B-cell","B-cell","B-cell","B-cell",
"AML","AML","AML",
"AML","AML","AML","AML","AML","AML","AML","AML")
cellcol<-rep("red",length(cells))
cellcol[cells=="B-cell"]="blue"
cellcol[cells=="T-cell"]="green"

#note remove endogenous controls- it does make a difference!


#Plots
golub1<-t(golub[-c(1:20),])
goluball<-GE.biplot(golub1,tit="GE-Biplot: all Golub Data",clabs=cells,cclr=cellcol,gclr="black",
opt.log=F,opt.stand=F,cex.lab=.8)
# golubchip<-GE.ChipPlot(golub1,tit="Chip-plot: all Golub Data", clabs=cells, gclr="red",cclr=cellcol,
# opt.log=F, opt.stand=F)
# tit1<-"Top 2.5% most varying Genes"
# lgmad<-apply(golub1, 2, mad)
# setmd<-lgmad>quantile(lgmad,.975)
# golubmad<-GE.biplot(golub1[,setmd],tit="GE-biplot: MAD selection (2.5%)",clabs=cells,cclr=cellcol,gclr="black",
# opt.log=F,opt.stand=F,cex.lab=.8)
# gm<-by(golub1, INDICES=list(cells), FUN=mean)
# genemeans<-cbind(gm[[2]],gm[[3]],gm[[1]]) # rearrange to order as  AML, B-cell and T-cell
# golubmad<-GE.GenePlot(golub1[,setmd], tit="Gene-plot : MAD selection (2.5%)",icon="GE.Icon", genedata=genemeans[setmd,],gclr="black",
# axes=c(1,2), opt.gf=T, opt.log=F, opt.stand=F,cex.lab=.8)
#GE.identifyGenes(golubmad$Genes,(golub.gnames[-c(1:20),2][setmd]))
}

#us age
par(mfrow=c(2,2),mar=c(2,2,2,2),ask=T)
GE.demo()
