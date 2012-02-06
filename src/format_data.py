'''
Format the data in a format that corresponds with the algorithm we
will be implementing
'''

'''
"Let X be the preprocessed gene expression data matrix, with n microarray samples
in the rows and p variable genes in the columns."

Our data should be stored as:

        0        1        2        ...        p
        
0       x_{0,0}  x_{0,1}  x_{0,2}             x_{0,p}
1       x_{1,0}  x_{1,1}  x_{1,2}             x_{1,p}  
...     
n       x_{n,0}  x_{n,1}  x_{n,2}             x_{n,p}


The gene number will be unique labeled 0 to p. A corresponding data structure
will be created to look up what gene the columns represent.
'''

'''
Matrix G is an n x r matrix representing microarray samples
Matrix H is an p x r matrix representing variable gene effects
'''

import os
import csv
import numpy as np

UNFORMATTED_PATH  = "../data/colon tumor/raw unformatted data/"
FORMATTED_PATH    = "../data/colon tumor/raw formatted data/"
UNFORMATTED_INPUT = "I2000.txt"
FORMATTED_OUTPUT  = "I2000.csv"
SAMPLES_COUNT     = 62
GENES_COUNT       = 2000

# Initialize list of lists
samples_list = [[]  for i in range(SAMPLES_COUNT)]



def check_paths():
    unformatted_path_exists = os.path.exists(UNFORMATTED_PATH)
    assert unformatted_path_exists, "Error: " + UNFORMATTED_PATH + " does not exist!"

    formatted_path_exists = os.path.exists(FORMATTED_PATH)
    assert formatted_path_exists, "Error: " + FORMATTED_PATH + " does not exist!"


def read_raw_I2000():
    I2000_file = open(UNFORMATTED_PATH + UNFORMATTED_INPUT,'r')
    
    # Read in the I2000 gene expression samples
    while True:
        line = I2000_file.readline()
        
        if not line:
            break
        
        # Handle whitespace
        if not line.strip():
            continue
        
        single_gene_col = line.split()
        
        assert len(single_gene_col) == 62, "Error: Data should have 62 samples per gene entry!" 
        
        for s in range(len(single_gene_col)):
            samples_list[s].append(np.double(single_gene_col[s]))
            
    # Verify data was formatted as expected, namely, 2000 genes per sample
    for s in range(SAMPLES_COUNT):
        assert len(samples_list[s]) == GENES_COUNT, "Error: Sample " + s + " has " + len(samples_list[s]) + " genes instead of " + GENES_COUNT

    print "Success: I2000 file read in."


def write_formatted_I2000():
    # Note: Write in binary mode so no extraneous new lines appear
    writer =  csv.writer(open(FORMATTED_PATH + FORMATTED_OUTPUT, 'wb'))
    
    for s in range(len(samples_list)):
        writer.writerow(samples_list[s])

    print "Success: I2000 comma seperated value file written"


read_raw_I2000()
write_formatted_I2000()






    