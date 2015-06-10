Pipeline containing a set of tools processing Next Generation Sequence data obtained from Degenerate Oligonucleotide Primed PCR libraries of isolated (flow sorted or microdissected) chromosomes to get some statistics and graphs. Pipeline has been designed for B chromosome data analysis but may be used for any sequence data. Starting point and input for this pipeline is paired end fastq files containing reads.

To start pipeline go to pipeline folder and use files located there. Additional information is in config file.

Dependancies:
1. cutadapt (tested on v.1.6)
2. bowtie2 (tested on v.2.1.0, 2.2.4)
3. Pysam python package (tested on v.0.8.1)
4. bedtools (tested on v.2.17.0)
5. DNAcopy R package ("http://bioconductor.org/biocLite.R")

Anolis directory holds modified and extended scripts of this pipeline used to analyze anolis sorted micro and sex chromosomes.
