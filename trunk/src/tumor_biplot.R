# Change to your respective trunk
TRUNK_DIR<-"C:/Users/King Charles/workspace/cap-6545-microarray-analysis/"
# Relative to your trunk directory
I2000_PATH<-"data/colon tumor/raw formatted data/I2000.csv"

# Sigh, R does not have good support for telling the path of the file that you're executing
setwd(TRUNK_DIR)
# Use the original GE-Biplot library until you've got a good handle on it.
source("public/ge-biplot.r", chdir=TRUE)

# Read in formatted tumor data
tumor_data <- read.csv(file=I2000_PATH,header=FALSE,sep=",");
# We need to do a matrix transpose to get it into the format it wants:
# rows = genes, cols = samples
tumor_data <- t(tumor_data)

# TODO: Read in tissue samples from formatted directory