# inScan

inScan is devoloping for finding genomoic insertion variation from long 
reads(from both Pacbio and Nanopore sequencing).
inScan differs from the published stat of art Structure Variation detecting 
tools(Sniffles and NanoSV) for long reads sequencing technology in two ways. First, 
inScan can find complex insertions when the insert sequence mapped to another 
chromosome. Second, inScan detects insertions in a given region, therefore useful 
for quickly checking if there are insertions in the region(Short Tandem Repeat region
for example) of interest.


# Getting started

```sh
git clone https://github.com/Archieyoung/inScan.git
cd inScan
python3 inScan.py input.bam input.bed output.json
```

