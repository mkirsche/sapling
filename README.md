SAPLING: Suffix Array Piecewise Linear INdex for Genomics  
  
A method for achieving faster suffix array queries compared to binary search.  
Currently the tool is still in development.  
  
Building:  
./build.sh  
  
Running SAPLING:  
./sapling <genome>  <suffix array file>  <sapling file>  
    genome is a fasta file containing the genome to query  
    suffix array file is where the suffix array will be written to  
    sapling file is where the sapling data structure will be written to  
  
Running binary search:  
./binarysearch <genome>  <suffix array file>  
    genome is a fasta file containing the genome to query  
    suffix array file is where the suffix array will be written to 
