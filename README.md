# WaIS
Where are Insertion Sequences

A package to find where the insertion sequences (IS) are located in a reference genome, or an assembly - using short-read sequences. Additional to finding the point of insertion, the orientation of the insertion and the IS type is also reported.


## Installation

[Installation instructions are here.](Installation.md)


## WaIS

Required parameters: 
```sh 
--outputDir
--ISseqs
--reads_1
--reads_2 
```

And, one of the following:
```sh
--runSpades
--assembly
```

### Running 
```sh

# Invoke help 
python3 wais.py --help 

# When you have already assembled the reads (we recommend using SPAdes, it takes longer than SKEASA, but provides higher accuracy for identifying IS insertions).
python3 wais/wais.py --outputDir SRR9988840 --ISseqs Example/IS_seqs_found_in_BP.fasta --reads_1 Example/SRR9988840_1.fastq.gz --reads_2 Example/SRR9988840_2.fastq.gz --assembly Example/contigs.fasta 

# Alternatively, also run spades with WaIS. 
python3 ... 

# Map the identified insertions to a reference genome.
python3 ... 

```

### Outputs

