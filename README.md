# SAPLING: Suffix Array Piecewise Linear INdex for Genomics  
  
A method for achieving faster suffix array queries compared to binary search.  
  
### Building:
```
cd src
./build.sh
make
```
  
### Running SAPLING test:  
```
./sapling_example <genome>  [<suffix array file>  <sapling file>]
    genome is a fasta file containing the genome to query  
    suffix array file is where the suffix array will be written to  
    sapling file is where the sapling data structure will be written to  
```
  
### Running aligner:
```
./align <query> <ref> <outfile>
    query is a FASTQ file containing reads
    ref is a FASTA file containing the reference genome
    outfile is where to write the SAM-format alignments
```
    
  
### Running binary search: 
```
./binarysearch <genome>  <suffix array file>  
    genome is a fasta file containing the genome to query  
    suffix array file is where the suffix array will be written to 
```
