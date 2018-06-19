# inScan

inScan is devoloping for finding genomoic insertion variation from long 
reads(from both Pacbio and Nanopore sequencing).
inScan differs from the published state of art Structure Variation detecting 
tools(Sniffles and NanoSV) for long reads sequencing technology in two ways. First, 
inScan can find complex insertions when the insert sequence mapped to another 
chromosome. Second, inScan detects insertions in a given region, therefore useful 
for quickly checking if there are insertions in the region(Short Tandem Repeat region
for example) of interest.

# Dependency

Python3.6.2 or later

* pysam

# Getting started

```sh
git clone git@github.com:Nextomics/InScan.git
cd inScan
python3 inScan.py input.bam input.bed output.json
```

inScan takes three positional arguments, a bam file, a bed file, a output file name.

inScan has been tested using NGMLR, BWA mem, Minimap2 output bam files. 
Generally, inScan will work for a bam file with "SA" tag.

bed file contains the regions to be tested.

output json: {"region":{"reads":[["chromosome","start","end","insert_size"],...],...},...}







