A pipeline for BCR repertoire libraries from  - UMI Barcoded Illumina MiSeq 325+275 paired-end.


Library preperation and sequencing method:

The sequences were amplified specific primers 1.constant region human Primers 2. V region primers. 
The generated libraries were then sequenced with Illumins MiSeq 325+275.

Input files:

* Pair-end reads Sample_R1.fastq and Sample_R2.fastq 
* primers sequences - constant region specific , v specific primers.
* Assemble pairs reference
* Isotype specific primers

Output file:

1. Sample.fasta - Two processed fasta files, one for Heavy chain sequence, and one for Light
2. log tab file for each steps
3. report for some of the steps


Pipeline container:

* Docker: immcantation/suite:4.3.0


Sequence processing steps:

1. FilterSeq quality
2. MaskPrimer align
3. PairAwk
4. AlignSeqts muscle
5. BuildConsensus
6. PairAwk
7. AssemblePairs de novo
8. AssemblePairs reference
9. MaskPrimer isotype
10. FilterSeq quality
11. FilterSeq length
12. CollapseSeq
13. SplitSeq to heavy and light

