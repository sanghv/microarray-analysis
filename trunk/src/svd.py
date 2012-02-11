'''
Created on Jan 29, 2012

@author: King Charles
edit: cw 

'''
from scipy import linalg, mat;

# code for svd
def svd():
	A = mat( [[1,0,0,0,2],
          	[0,0,3,0,0],
          	[0,0,0,0,0],
          	[0,4,0,0,0]] )

	print "================="
	print "A:\n",  A
	print "=================\n"

	# S for sigma
	U, S, V = linalg.svd( A )

	# This one does not match the worked out wikipedia tutorial. Hmm.
	print "U is a matrix whose columns are the eigenvectors of the A*A^T matrix."
	print "These are termed the left eigenvectors.\n"
	print "U:\n", U, "\n"

	print "S is a matrix whose diagonal elements are the singular values of A."
	print "This is a diagonal matrix, so its nondiagonal elements are zero by definition.\n"
	print "Sigma:\n", S, "\n"

	print "V is a matrix whose columns are the eigenvectors of the ATA matrix."
	print "These are termed the right eigenvectors. VT is the transpose of V.\n"
	print "V^T:\n", V, "\n"

# main function used for test functions/classes
def main():
	svd()	

# main used for test cases: pust test calls in main()
if __name__ == "__main__":
	main()
