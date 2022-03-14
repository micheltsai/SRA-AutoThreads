# Benga
Bacterial Epidemiology NGs Analysis (BENGA) framework and pipeline.

# Requirements
* python 3.6+
  * biopython
  * scipy
  * numba
  * numpy
  * pandas
  * fastcluster
  * matplotlib
* ncbi-blast
* diamond
* prodigal
* prokka
* roary

# Usage
***cgMLST profiling***
```
profiling -i fasta_file -o result.tsv --scheme scheme.faa --prodigaltf train_file.trn
```
***cgMLST makedb***
```
makedb.py -i input_path -o output_path
```

***draw dendrogram***
```
dendrogram.py -i profile_1 profile_2 profile_3 ... -o output_path
```
