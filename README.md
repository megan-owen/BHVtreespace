This repository contains the latest source code and jar files for various computations in the Billera-Holmes-Vogtmann treespace, including distance, mean, and variance.

See [http://comet.lehman.cuny.edu/owen/code.html](http://comet.lehman.cuny.edu/owen/code.html) for links to corresponding journal papers and manuels for using the code.

## Addendums to manuels
### GTP
March 21, 2023: Added an option (`-f <otherTreeFile>`) for GTP to read in two files and compute either all BHV distances between the trees in each file (default) or to compute only the BHV distance betweeen the i-th trees in each file (also use `-p` option).  

Examples:
1. `java -jar -u -f treeFile2.txt -o outFile.txt treeFile1.txt`
  Computes BHV distance between all unrooted trees in treeFile1.txt and all unrooted trees in treeFile2.txt, storing the distances in the file outFile.txt.
  
2. `java -jar -f treeFile2.txt -p -o outFile.txt treeFile1.txt`
Computes BHV distance between the first tree in treeFile1.txt and the first tree in treeFile2.txt, the second tree in treeFile1.txt and the second tree in treeFile2.txt, etc.  Stores the distances in file outFile.txt and assumes all trees are rooted.



