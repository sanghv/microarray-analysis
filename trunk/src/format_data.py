import os
import csv
import numpy as np
from operator import itemgetter

'''
Dataset
This is the base class we will inherit from and implement reading and writing
the microarray data into a standard format that we will later use in our
KPCA-Biplot application.
'''
class DataSet:

    name = "Dataset" # Override this

    expression_profiles = None
    classification_profiles = None
    gene_list = None

    def verify(self):
        pass

    def read(self):
        pass

    def write(self):
        pass


class TumorDataSet(DataSet):

    '''
    Data set constants
    '''
    UNFORMATTED_PATH         = "../data/COLON/unformatted/"
    FORMATTED_PATH           = "../data/COLON/formatted/"

    UNFORMATTED_I2000_INPUT  = "I2000.txt" # contains gene expression profiles
    UNFORMATTED_TISSUE_INPUT = "tissues.txt" # contains classifications
    UNFORMATTED_GENES_INPUT  = "names.txt" # contains gene names

    FORMATTED_PROFILE_OUTPUT = "expression_profiles.csv"  # contains genes and samples
    FORMATTED_CLASSIFICATION_OUTPUT = "classification.txt" # contains list of classifications
    FORMATTED_GENE_OUTPUT = "genes.txt" # contains list of genes

    SAMPLES_COUNT     = 62
    GENES_COUNT       = 2000

    name = "Tumor Data Set"

    def __init__(self):
        # Initialize list of lists to be used for reading and writing
        self.expression_profiles = [[]  for i in range(self.SAMPLES_COUNT)]
        self.classification_profiles  = []
        self.gene_list = []

    def verify(self):
        unformatted_path_exists = os.path.exists(self.UNFORMATTED_PATH)
        assert unformatted_path_exists, "Error: " + self.UNFORMATTED_PATH + " does not exist!"

        formatted_path_exists = os.path.exists(self.FORMATTED_PATH)
        assert formatted_path_exists, "Error: " + self.FORMATTED_PATH + " does not exist!"

        return unformatted_path_exists and formatted_path_exists

    def read(self):
        I2000_read_was_success  = self.read_I2000()
        tissue_read_was_success = self.read_tissues()
        gene_read_was_success   = self.read_genes()

        return I2000_read_was_success and tissue_read_was_success and gene_read_was_success

    def read_I2000(self):
        I2000_file = open(self.UNFORMATTED_PATH + self.UNFORMATTED_I2000_INPUT,'r')

        # Read in the I2000 gene expression samples
        while True:
            line = I2000_file.readline()

            if not line:
                break

            # Handle whitespace
            if not line.strip():
                continue

            single_gene_col = line.split()

            #TODO DON'T USE ACTUAL NUMBERS, USE VARIABLES
            correct_sample_len = len(single_gene_col) == self.SAMPLES_COUNT
            assert correct_sample_len, "Error: Data should have " + str(self.SAMPLES_COUNT) + \
                                       " samples per gene entry!"

            for s in range(len(single_gene_col)):
                self.expression_profiles[s].append(np.double(single_gene_col[s]))

        # Verify data was formatted as expected, namely, 2000 genes per sample
        for s in range(self.SAMPLES_COUNT):
            correct_gene_len = len(self.expression_profiles[s]) == self.GENES_COUNT
            assert correct_gene_len, "Error: Sample " + str(s) + " has " + str(len(self.expression_profiles[s])) + \
                                     " genes instead of " + str(self.GENES_COUNT)

        print "Success: I2000 file read in."
        return True

    def read_tissues(self):
        tissue_file = open(self.UNFORMATTED_PATH + self.UNFORMATTED_TISSUE_INPUT,'r')

        # Read in the tissue classifications
        while True:
            tissue_sample = tissue_file.readline()

            if not tissue_sample:
                break

            # Handle whitespace
            if not tissue_sample.strip():
                continue

            # TODO: any error checking necessary?
            tissue_value = np.int(tissue_sample)

            if tissue_value >= 0:
                self.classification_profiles.append('+') # Sample came from a patient with a colon tumor
            else:
                self.classification_profiles.append('-') # Sample came from a patient with a healthy colon

        print "Success: tissue file read in."
        return True

    def read_genes(self):
        gene_file = open(self.UNFORMATTED_PATH + self.UNFORMATTED_GENES_INPUT,'r')

        # Read in the tissue classifications
        while True:
            gene = gene_file.readline()

            if not gene:
                break

            # Handle whitespace
            if not gene.strip():
                continue

            gene_info = gene.split()

            GENE_ALT_INDEX = 0
            GENE_INDEX     = 1
            assert len(gene_info) > 1, "Error: not enough fields in gene_info"
            gene_name = gene_info[GENE_INDEX]

            if gene_name.upper() == "CONTROL":
                gene_name = gene_info[GENE_ALT_INDEX]
                self.gene_list.append(gene_name)
            else:
                assert len(gene_name) > 1, "Error: not enough chars in gene name"
                assert gene_name[0].isalpha(), "Error: not in the correct format, " + gene_name
                assert gene_name[1:].isdigit(), "Error: not in the correct format, " + gene_name

                self.gene_list.append(gene_name)


        print "Success: gene file read in."
        return True


    def write(self):
        expression_write_was_success     = self.write_expression_profile()
        classification_write_was_success = self.write_classifications()
        gene_write_was_success           = self.write_genes()

        return expression_write_was_success and classification_write_was_success and gene_write_was_success

    def write_expression_profile(self):
        assert len(self.expression_profiles) == self.SAMPLES_COUNT, "Error: profile count differs from expected"
        # Note: Write in binary mode so no extraneous new lines appear
        expression_writer =  csv.writer(open(self.FORMATTED_PATH + self.FORMATTED_PROFILE_OUTPUT, 'wb'))

        for s in range(len(self.expression_profiles)):
            expression_writer.writerow(self.expression_profiles[s])

        print "Success: gene expression file written"
        return True

    def write_classifications(self):
        assert len(self.classification_profiles) == self.SAMPLES_COUNT, "Error: Sample count differs from expected"
        classification_file = open(self.FORMATTED_PATH + self.FORMATTED_CLASSIFICATION_OUTPUT, 'w')

        # Write each tissue sample on it's own line
        for classification_sample in self.classification_profiles:
            classification_file.write(classification_sample + "\n")

        print "Success: classification file written"
        return True

    def write_genes(self):
        assert len(self.gene_list) == self.GENES_COUNT, "Error: Gene count differs from expected"
        gene_file = open(self.FORMATTED_PATH + self.FORMATTED_GENE_OUTPUT, 'w')

        # Write each tissue sample on it's own line
        for gene in self.gene_list:
            gene_file.write(gene + "\n")

        print "Success: gene file written"
        return True


class LymphomaDataSet(DataSet):

    name = "Lymphoma Data Set"

    UNFORMATTED_PATH         = "../data/LYMPHOMA/unformatted/"
    FORMATTED_PATH           = "../data/LYMPHOMA/formatted/"

    UNFORMATTED_DATASET_INPUT  = "dataset.txt" # contains genes, samples, and classification
    FORMATTED_PROFILE_OUTPUT = "expression_profiles.csv"  # contains genes and samples
    FORMATTED_CLASSIFICATION_OUTPUT = "classification.txt" # contains list of classifications
    FORMATTED_GENE_OUTPUT = "genes.txt" # contains list of genes

    SAMPLES_COUNT = 77
    GENES_COUNT   = 2647 # Page says 7129 though, even though there data set only has 2647

    def __init__(self):
        # Initialize list of lists to be used for reading and writing
        self.expression_profiles = [[]  for i in range(self.SAMPLES_COUNT)]
        self.classification_profiles = []
        self.gene_list = []

    def verify(self):
        unformatted_path_exists = os.path.exists(self.UNFORMATTED_PATH)
        assert unformatted_path_exists, "Error: " + self.UNFORMATTED_PATH + " does not exist!"

        formatted_path_exists = os.path.exists(self.FORMATTED_PATH)
        assert formatted_path_exists, "Error: " + self.FORMATTED_PATH + " does not exist!"

        return unformatted_path_exists and formatted_path_exists

    def read(self):
        dataset_file = open(self.UNFORMATTED_PATH + self.UNFORMATTED_DATASET_INPUT,'r')

        # Read in the gene expression profiles
        # NOTE: Format is:
        # <Gene Name_{1}> <Samples_{1 .. SAMPLES_COUNT}>
        # ...
        # <Gene Name_{GENES_COUNT}> <Samples_{1 .. SAMPLES_COUNT}>
        # <ydlbcl> <CLASSIFICATION_{1 .. SAMPLES_COUNT} where 0 = {?} and 1 = {?}>
        while True:
            line = dataset_file.readline()

            if not line:
                break

            # Handle whitespace
            if not line.strip():
                continue

            # Read classifications
            if line.startswith("ydlbcl"):
                classifications = line.split()

                correct_classifications_len = len(classifications) == self.SAMPLES_COUNT + 1
                assert correct_classifications_len, "Error: Data should have " + str(self.SAMPLES_COUNT) + \
                                                    " classifications!"
                # remove identifier "ydlbcl" from the list
                del classifications[0]

                self.classification_profiles += classifications

                continue
            # Else, read gene intensities for single gene
            gene_intensity_samples = line.split()

            correct_sample_len = len(gene_intensity_samples) == self.SAMPLES_COUNT + 1
            assert correct_sample_len, "Error: Data should have " + str(self.SAMPLES_COUNT) + \
                                        " samples per gene entry!"

            # remove gene name from the list
            gene_name = gene_intensity_samples[0]
            del gene_intensity_samples[0]
            self.gene_list.append(gene_name)

            for s in range(len(gene_intensity_samples)):
                self.expression_profiles[s].append(np.double(gene_intensity_samples[s]))

        # Verify data was formatted as expected, namely, <GENES_COUNT> genes per sample
        for s in range(self.SAMPLES_COUNT):
            correct_gene_len = len(self.expression_profiles[s]) == self.GENES_COUNT
            assert correct_gene_len, "Error: Sample " + str(s) + " has " + str(len(self.expression_profiles[s])) + \
                                     " genes instead of " + str(self.GENES_COUNT)

        # Verify gene list has correct amount of genes
        correct_gene_list_len = len(self.gene_list) == self.GENES_COUNT
        assert correct_gene_list_len, "Error: Gene list has" + str(len(self.gene_list)) + " genes" + \
                                      "instead of " + str(self.GENES_COUNT)

        # Verify classification has correct amount of samples
        correct_classifications_len = len(self.classification_profiles) == self.GENES_COUNT
        assert correct_gene_list_len, "Error: Gene list has" + str(len(self.classification_profiles)) + \
                                      " classifications" +"instead of " + str(self.GENES_COUNT)

        print "Success: dataset, gene names, and classifications read in"
        return True



    def write(self):

        # Note: Write in binary mode so no extraneous new lines appear
        gene_expression_writer =  csv.writer(open(self.FORMATTED_PATH + self.FORMATTED_PROFILE_OUTPUT, 'wb'))

        for s in range(len(self.expression_profiles)):
            gene_expression_writer.writerow(self.expression_profiles[s])


        gene_file = open(self.FORMATTED_PATH + self.FORMATTED_GENE_OUTPUT, 'w')

        # Write each gene name on its own line
        for gene in self.gene_list:
            gene_file.write(gene + "\n")


        classification_file = open(self.FORMATTED_PATH + self.FORMATTED_CLASSIFICATION_OUTPUT, 'w')

        # Write each classification sample on its own line
        for classification_sample in self.classification_profiles:
            classification_file.write(classification_sample + "\n")

        print "Success: Gene expression profile, classifications, and gene names written"
        return True


class ProstateDataSet(DataSet):

    name = "Prostate Data Set"

    UNFORMATTED_PATH         = "../data/PROSTATE/unformatted/"
    FORMATTED_PATH           = "../data/PROSTATE/formatted/"

    UNFORMATTED_DATASET_INPUT  = "dataset.txt" # contains genes, samples, and classification
    FORMATTED_PROFILE_OUTPUT = "expression_profiles.csv"  # contains genes and samples
    FORMATTED_CLASSIFICATION_OUTPUT = "classification.txt" # contains list of classifications
    FORMATTED_GENE_OUTPUT = "genes.txt" # contains list of genes

    SAMPLES_COUNT = 102
    GENES_COUNT   = 2135

    def __init__(self):
        # Initialize list of lists to be used for reading and writing
        self.expression_profiles = [[]  for i in range(self.SAMPLES_COUNT)]
        self.classification_profiles = []
        self.gene_list = []

    def verify(self):
        unformatted_path_exists = os.path.exists(self.UNFORMATTED_PATH)
        assert unformatted_path_exists, "Error: " + self.UNFORMATTED_PATH + " does not exist!"

        formatted_path_exists = os.path.exists(self.FORMATTED_PATH)
        assert formatted_path_exists, "Error: " + self.FORMATTED_PATH + " does not exist!"

        return unformatted_path_exists and formatted_path_exists

    def read(self):
        dataset_file = open(self.UNFORMATTED_PATH + self.UNFORMATTED_DATASET_INPUT,'r')

        # Read in the gene expression profiles
        # NOTE: Format is:
        # <Gene Name_{1}> <Samples_{1 .. SAMPLES_COUNT}>
        # ...
        # <Gene Name_{GENES_COUNT}> <Samples_{1 .. SAMPLES_COUNT}>
        # <ysingh> <CLASSIFICATION_{1 .. SAMPLES_COUNT} where 0 = {?} and 1 = {?}>
        while True:
            line = dataset_file.readline()

            if not line:
                break

            # Handle whitespace
            if not line.strip():
                continue

            # Read classifications
            if line.startswith("ysingh"):
                classifications = line.split()

                correct_classifications_len = len(classifications) == self.SAMPLES_COUNT + 1
                assert correct_classifications_len, "Error: Data should have " + str(self.SAMPLES_COUNT) +\
                                                    " classifications!"
                # remove identifier "ydlbcl" from the list
                del classifications[0]

                self.classification_profiles += classifications

                continue
                # Else, read gene intensities for single gene
            gene_intensity_samples = line.split()

            correct_sample_len = len(gene_intensity_samples) == self.SAMPLES_COUNT + 1
            assert correct_sample_len, "Error: Data should have " + str(self.SAMPLES_COUNT) +\
                                       " samples per gene entry!"

            # remove gene name from the list
            gene_name = gene_intensity_samples[0]
            del gene_intensity_samples[0]
            self.gene_list.append(gene_name)

            for s in range(len(gene_intensity_samples)):
                self.expression_profiles[s].append(np.double(gene_intensity_samples[s]))

        # Verify data was formatted as expected, namely, <GENES_COUNT> genes per sample
        for s in range(self.SAMPLES_COUNT):
            correct_gene_len = len(self.expression_profiles[s]) == self.GENES_COUNT
            assert correct_gene_len, "Error: Sample " + str(s) + " has " + str(len(self.expression_profiles[s])) +\
                                     " genes instead of " + str(self.GENES_COUNT)

        # Verify gene list has correct amount of genes
        correct_gene_list_len = len(self.gene_list) == self.GENES_COUNT
        assert correct_gene_list_len, "Error: Gene list has" + str(len(self.gene_list)) + " genes" +\
                                      "instead of " + str(self.GENES_COUNT)

        # Verify classification has correct amount of samples
        correct_classifications_len = len(self.classification_profiles) == self.GENES_COUNT
        assert correct_gene_list_len, "Error: Gene list has" + str(len(self.classification_profiles)) +\
                                      " classifications" +"instead of " + str(self.GENES_COUNT)

        print "Success: dataset, gene names, and classifications read in"
        return True

    def write(self):
        # Note: Write in binary mode so no extraneous new lines appear
        gene_expression_writer =  csv.writer(open(self.FORMATTED_PATH + self.FORMATTED_PROFILE_OUTPUT, 'wb'))

        for s in range(len(self.expression_profiles)):
            gene_expression_writer.writerow(self.expression_profiles[s])


        gene_file = open(self.FORMATTED_PATH + self.FORMATTED_GENE_OUTPUT, 'w')

        # Write each gene name on its own line
        for gene in self.gene_list:
            gene_file.write(gene + "\n")


        classification_file = open(self.FORMATTED_PATH + self.FORMATTED_CLASSIFICATION_OUTPUT, 'w')

        # Write each classification sample on its own line
        for classification_sample in self.classification_profiles:
            classification_file.write(classification_sample + "\n")

        print "Success: Gene expression profile, classifications, and gene names written"
        return True


class LeukemiaDataSet(DataSet):

    '''
    Data set constants
    '''
    UNFORMATTED_PATH         = "../data/LEUKEMIA/unformatted/"
    FORMATTED_PATH           = "../data/LEUKEMIA/formatted/"

    UNFORMATTED_TEST_DATA_INPUT           = "golub.test"
    UNFORMATTED_TEST_CLASSIFICATION_INPUT = "golub.test.pheno"

    UNFORMATTED_TRAIN_DATA_INPUT           = "golub.train"
    UNFORMATTED_TRAIN_CLASSIFICATION_INPUT = "golub.train.pheno"

    FORMATTED_PROFILE_OUTPUT = "expression_profiles.csv"  # contains genes and samples
    FORMATTED_CLASSIFICATION_OUTPUT = "classification.txt" # contains list of classifications
    FORMATTED_GENE_OUTPUT = "genes.txt" # contains list of genes

    SAMPLES_TRAIN_COUNT = 38
    SAMPLES_TEST_COUNT  = 34
    SAMPLES_TOTAL_COUNT = SAMPLES_TRAIN_COUNT + SAMPLES_TEST_COUNT
    GENES_COUNT       = 7129

    name = "Leukemia Data Set"
    expression_train_profiles = None
    expression_test_profiles  = None
    verification_gene_list    = None
    test_classifications      = None
    train_classifications     = None

    def __init__(self):
        # Initialize list of lists to be used for reading and writing

        #self.expression_profiles = [[]  for i in range(self.SAMPLES_TOTAL_COUNT)]
        self.expression_train_profiles = [[]  for i in range(self.SAMPLES_TRAIN_COUNT)]
        self.expression_test_profiles  = [[]  for i in range(self.SAMPLES_TEST_COUNT)]

        self.test_classifications     = {} # dictionary because the labeled samples are out of order
        self.train_classifications    = {} # dictionary because the labeled samples are out of order
        self.classification_profiles  = [] # they will be sorted into this list

        self.gene_list = []
        self.verification_gene_list = [] # Since we have to read the gene list twice, why not verify they're the same


    def verify(self):
        unformatted_path_exists = os.path.exists(self.UNFORMATTED_PATH)
        assert unformatted_path_exists, "Error: " + self.UNFORMATTED_PATH + " does not exist!"

        formatted_path_exists = os.path.exists(self.FORMATTED_PATH)
        assert formatted_path_exists, "Error: " + self.FORMATTED_PATH + " does not exist!"

        return unformatted_path_exists and formatted_path_exists

    def read(self):
        train_read_was_success = self.read_train()
        test_read_was_success  = self.read_test()
        train_class_read_was_success = self.read_train_classifications()
        test_class_read_was_success  = self.read_test_classifications()

        return train_read_was_success and test_read_was_success and \
               train_class_read_was_success and test_class_read_was_success

    def read_train_classifications(self):
        dataset_file = open(self.UNFORMATTED_PATH + self.UNFORMATTED_TRAIN_CLASSIFICATION_INPUT,'r')

        # Read in the gene expression profiles
        # NOTE: Format is:
        # <Sample_{index=start}>..<Sample_{index=end}>
        # <Sample Number> <Classification> <Other Stuff>
        # ...
        while True:
            line = dataset_file.readline()

            if not line:
                break

            # Handle whitespace
            if not line.strip():
                continue

            # Skip headers
            if line.startswith("Samples"):
                continue

            # Else, read gene intensities for single gene
            labeled_sample = line.split()

            assert len(labeled_sample) == 11, "Error: Sample fields differ from expected"

            sample_num = np.int(labeled_sample[0])
            sample_classification = labeled_sample[1]

            assert not self.train_classifications.has_key(sample_num), "Error: duplicate key values"

            self.train_classifications[sample_num] = sample_classification


        assert len(self.train_classifications) == self.SAMPLES_TRAIN_COUNT, "Error: Incorrect amount of samples read"

        print "Success: Train classifications read in successfully"
        return True

    def read_test_classifications(self):
        dataset_file = open(self.UNFORMATTED_PATH + self.UNFORMATTED_TEST_CLASSIFICATION_INPUT,'r')

        # Read in the gene expression profiles
        # NOTE: Format is:
        # <Sample_{index=start}>..<Sample_{index=end}>
        # <Sample Number> <Classification> <Other Stuff>
        # ...
        while True:
            line = dataset_file.readline()

            if not line:
                break

            # Handle whitespace
            if not line.strip():
                continue

            # Skip headers
            if line.startswith("Samples"):
                continue

            # Else, read gene intensities for single gene
            labeled_sample = line.split()

            assert len(labeled_sample) == 11, "Error: Sample fields differ from expected"

            sample_num = np.int(labeled_sample[0])
            sample_classification = labeled_sample[1]

            assert not self.test_classifications.has_key(sample_num), "Error: duplicate key values"

            self.test_classifications[sample_num] = sample_classification


        assert len(self.test_classifications) == self.SAMPLES_TEST_COUNT, "Error: Incorrect amount of samples read"

        print "Success: Test classifications read in successfully"
        return True


    '''
    THIS IS READ IN FIRST
    '''
    def read_train(self):
        dataset_file = open(self.UNFORMATTED_PATH + self.UNFORMATTED_TRAIN_DATA_INPUT,'r')

        # TODO UPDATE COMMENT V V V
        # Read in the gene expression profiles
        # NOTE: Format is:
        # <Sample_{index=start}>..<Sample_{index=end}>
        # <Gene Name_{1}> <Samples_{1 .. SAMPLES_COUNT}>
        # ...
        # <Gene Name_{GENES_COUNT}> <Samples_{1 .. SAMPLES_COUNT}>
        while True:
            line = dataset_file.readline()

            if not line:
                break

            # Handle whitespace
            if not line.strip():
                continue

            # Skip sample number indices line since they're already in order
            if line.startswith("X1 X2"):
                continue

            # Else, read gene intensities for single gene
            gene_intensity_samples = line.split()

            correct_sample_len = len(gene_intensity_samples) == self.SAMPLES_TRAIN_COUNT + 1
            assert correct_sample_len, "Error: Train Data should have " + str(self.SAMPLES_TRAIN_COUNT) +\
                                       " samples per gene entry!"

            # remove gene name from the list
            gene_name = gene_intensity_samples[0]
            del gene_intensity_samples[0]

            if len(self.gene_list) < self.GENES_COUNT:
                self.gene_list.append(gene_name)

            for s in range(len(gene_intensity_samples)):
                self.expression_train_profiles[s].append(np.double(gene_intensity_samples[s]))

        # Verify data was formatted as expected, namely, <GENES_COUNT> genes per sample
        for s in range(self.SAMPLES_TRAIN_COUNT):
            correct_gene_len = len(self.expression_train_profiles[s]) == self.GENES_COUNT
            assert correct_gene_len, "Error: Sample " + str(s) + " has " + \
                                     str(len(self.expression_train_profiles[s])) + \
                                     " genes instead of " + str(self.GENES_COUNT)

        # Verify gene list has correct amount of genes
        correct_gene_list_len = len(self.gene_list) == self.GENES_COUNT
        assert correct_gene_list_len, "Error: Gene list has" + str(len(self.gene_list)) + " genes" +\
                                      "instead of " + str(self.GENES_COUNT)

        print "Success: Train data and genes successfully read in"
        return True



    '''
    This is read in SECOND
    '''
    def read_test(self):
        dataset_file = open(self.UNFORMATTED_PATH + self.UNFORMATTED_TEST_DATA_INPUT,'r')

        # TODO UPDATE COMMENT V V V
        # Read in the gene expression profiles
        # NOTE: Format is:
        # <Gene Name_{1}> <Samples_{1 .. SAMPLES_COUNT}>
        # ...
        # <Gene Name_{GENES_COUNT}> <Samples_{1 .. SAMPLES_COUNT}>
        while True:
            line = dataset_file.readline()

            if not line:
                break

            # Handle whitespace
            if not line.strip():
                continue

            # Skip sample number indices line since they're already in order
            if line.startswith("X1 X2"):
                continue

            # Else, read gene intensities for single gene
            gene_intensity_samples = line.split()

            correct_sample_len = len(gene_intensity_samples) == self.SAMPLES_TEST_COUNT + 1
            assert correct_sample_len, "Error: Test Data should have " + str(self.SAMPLES_TEST_COUNT) +\
                                       " samples per gene entry!"

            # remove gene name from the list
            gene_name = gene_intensity_samples[0]
            del gene_intensity_samples[0]

            if len(self.verification_gene_list) < self.GENES_COUNT:
                self.verification_gene_list.append(gene_name)

            for s in range(len(gene_intensity_samples)):
                self.expression_test_profiles[s].append(np.double(gene_intensity_samples[s]))

        # Verify data was formatted as expected, namely, <GENES_COUNT> genes per sample
        for s in range(self.SAMPLES_TEST_COUNT):
            correct_gene_len = len(self.expression_test_profiles[s]) == self.GENES_COUNT
            assert correct_gene_len, "Error: Sample " + str(s) + " has " +\
                                     str(len(self.expression_test_profiles[s])) +\
                                     " genes instead of " + str(self.GENES_COUNT)


        # Verify gene list and verification gene list are the same
        assert len(self.gene_list) == len(self.verification_gene_list), "Error: Gene lengths differ between " + \
                                                                        "test and training data set"
        assert self.gene_list == self.verification_gene_list, "Error: Test and train data set contain different " + \
                                                              "genes"

        print "Success: Test data and genes successfully verified"
        return True


    def write(self):
        # Note: Write in binary mode so no extraneous new lines appear
        gene_expression_writer =  csv.writer(open(self.FORMATTED_PATH + self.FORMATTED_PROFILE_OUTPUT, 'wb'))

        # NOTE: BE CONSISTENT! COMBINE TRAIN DATA THEN TEST DATA

        # Combine test and train profiles
        self.expression_profiles = self.expression_train_profiles + self.expression_test_profiles

        for s in range(len(self.expression_profiles)):
            gene_expression_writer.writerow(self.expression_profiles[s])


        classification_file = open(self.FORMATTED_PATH + self.FORMATTED_CLASSIFICATION_OUTPUT, 'w')

        classification_file = open(self.FORMATTED_PATH + self.FORMATTED_CLASSIFICATION_OUTPUT, 'w')
        self.classification_profiles = sorted(self.test_classifications.items() +  self.train_classifications.items(), \
                                              key=itemgetter(0))
        self.classification_profiles = [value for key,value in self.classification_profiles] # get values from tuples

        # Write each tissue sample on it's own line
        for classification_item in self.classification_profiles:
            classification_file.write(classification_item + "\n")

        print "Success! Gene classifications and expression profiles written"

        return True



if __name__ == '__main__':

    datasets = [TumorDataSet(), LymphomaDataSet(), ProstateDataSet(), LeukemiaDataSet()]

    for i in range(len(datasets)):
        print "Dataset: " + datasets[i].name + "\n"
        assert datasets[i].verify(), "Verification of dataset " + datasets[i].name + " failed"
        assert datasets[i].read(), "Read of dataset " + datasets[i].name + " failed"
        assert datasets[i].write(), "Write of dataset " + datasets[i].name + " failed"
        print ""