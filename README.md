Pipeline containing set of tools designed to process Next Generation Sequence data to get some statistics and graphs. Pipeline has been designed for B chromosome data analysis but may be used for any sequence data. Starting point and input for this pipeline is paired fastq files containing reads.
To start pipeline simply run main shell script with .txt file as argument. That .txt file should contain name of samples and of reference genomes for those samples and should have following format:
1_name.control name_ref
2_name2.control name_ref2
3_name3 name_ref
(e.g.
1_CFA12.control canFam3
2_BTA26.control bosTau7
3_CPYB bosTau7
4_NPP6 canFam3
)
Note that controls should have .control after their name and alphabetical order of names must be same as alphabetical order of fastq files used.
To run this pipeline apart from fastq files you also should put following files in folder with this pipeline:
1. Bowtie2 indexed refernce genomes (named same as input .txt file) and indexed hg19 genome
2. Text .genome files containing information about number and sizes of chromosomes in reference genomes (named same as input .txt file)
3. Folder named bed.reg.files containing .bed files with locations of regions
For example if you wanted to launch this pipeline using input .txt file same as in example above you would need to have:
1. bosTau7.1.bt2, bosTau7.2.bt2, bosTau7.3.bt2, bosTau7.4.bt2, bosTau7.fa, bosTau7.rev.1.bt2, bosTau7.rev.2.bt2 files and same files for canFam3 and hg19 (named exactly)
2. bosTau7.genome, canFam3.genome
3. Folder named bed.reg.files containing 3_CPYB.reg.bed, 4_NPP6.reg.bed (note that controls dont need region files)

