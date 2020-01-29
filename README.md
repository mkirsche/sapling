# SAPLING: Suffix Array Piecewise Linear INdex for Genomics  
  
A method for achieving faster suffix array queries compared to binary search.
  
### Building:
```
cd src
make
```

### Pre-building the Suffix Array:
```
suffixarray/refToSuffixArray.sh <genome file (Fasta format)>
 ```
  
### Running SAPLING test:  
```
src/sapling_example <genome file (Fasta format)> 
  [saFn=<suffix array file (format described in manual, default <genome file>.sa)>] 
  [sapFn=<sapling file (format described in manual, default <genome file>.sap>] 
  [nb=<log number of buckets - overrides maxMem below, default None>] 
  [maxMem=<max number of buckets will be (genome size)/val, default 10>] 
  [k=<k, default 21>] 
  [nq=<number of queries, default 5000000>] 
  [errFn=<errors file if outputting them, default None>]
```
  
### Running aligner:
```
src/align <query (Fastq format)> <ref (Fasta format)> <outfile (Sam format)> 
  [num_seeds=<number of seeds to use for exact matching, default 7>]
  [sapling_k=<size of k-mers to use when building Sapling, default 16>]
  [flanking_sequence=<amount of padding to include when aligning region around seeds, default 2>]
  [max_hits=<maximum number of matches to try to extend per seed, default 32>]
```
    
  
### Running binary search: 
```
src/binarysearch <genome> <suffix array file>
    genome is a fasta file containing the genome to query
    suffix array file is where the suffix array will be written to
```

### Sapling Manual

Information about the files in this repository can be found in the [manual](https://github.com/mkirsche/sapling/wiki/Sapling-Manual).
