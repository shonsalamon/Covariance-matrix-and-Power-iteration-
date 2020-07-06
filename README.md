The cov program receives two filenames as command-line arguments. The 1st refers to the input file, and the 2nd is the output file. 
The input file is a matrix, and the output file produced should contain its covariance matrix.




The power iteration program receives two command-line arguments. The program receives three or four arguments:
1. Input matrix – file name of a symmetric matrix to perform power iteration on.
2. Initial vector – file name of an initial vector b0 to use for power iteration (optional).
3. Output file – file name to output the produced eigenvector into.
4. Implementation – a flag determining what sparse matrix implementation to use, either "-list" or "-array".

The program receives an input matrix and an implementation choice and outputs an estimate of its eigenvector
with the highest eigenvalue, by storing and manipulating the input matrix as a sparse matrix using the chosen implementation.
