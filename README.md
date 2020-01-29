# SAPLING: Suffix Array Piecewise Linear INdex for Genomics  
  
A method for achieving faster suffix array queries compared to binary search.
  
### Building:
```
cd src
./build.sh
make
```

### Pre-building the Suffix Array:
```
suffixarray/refToSuffixArray.sh <genome file (Fasta)>
 ```
  
### Running SAPLING test:  
```
src/sapling_example <genome file (Fasta)> 
  [saFn=<suffix array file (format described in manual)>] 
  [sapFn=<sapling file (format described in manual>] 
  [nb=<log number of buckets>] 
  [maxMem=<max number of buckets will be (genome size)/val>] 
  [k=<k>] 
  [nq=<number of queries>] 
  [errFn=<errors file if outputting them>]
```
  
### Running aligner:
```
src/align <query (Fastq)> <ref (Fasta)> <outfile (Sam)> 
  [num_seeds=<number of seeds to use for exact matching>]
  [sapling_k=<size of k-mers to use when building Sapling>]
  [flanking_sequence=<amount of padding to include when aligning region around seeds>]
  [max_hits=<maximum number of matches to try to extend per seed>]
```
    
  
### Running binary search: 
```
src/binarysearch <genome> <suffix array file>
    genome is a fasta file containing the genome to query
    suffix array file is where the suffix array will be written to
```

### Sapling Manual

Information about the files in this repository can be found in the [manual](https://github.com/mkirsche/sapling/wiki/Sapling-Manual).
