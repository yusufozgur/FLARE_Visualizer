# FLARE_Visualizer
Visualizes FLARE local ancestry inference using tagore.

Requires a specific version of tagore, which is found at https://github.com/Bahex/tagore, to be installed. Make sure tagore is present in PATH.

To see arguments, run the following command: Rscript FLARE_Visualizer.R --help

Files are produced by FLARE are taken as input, they are expected to be in the following format

```
##fileformat=VCFv4.2
##filedate=20230101
##source=flare.20Oct22.2a6.jar
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AN1,Number=1,Type=Integer,Description="Ancestry of first haplotype">
##FORMAT=<ID=AN2,Number=1,Type=Integer,Description="Ancestry of second haplotype">
##ANCESTRY=<ANCESTRY_first=0,ANCESTRY_second=1,ANCESTRY_third=2>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Person1 Person2 ...
1	100	.	G	T	.	PASS	.	GT:AN1:AN2	0|0:0:0	1|1:2:2 ...
...
```
