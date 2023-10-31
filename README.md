# FLARE_Visualizer
Visualizes FLARE local ancestry inference using tagore.

Requires a specific version of tagore, which is found at https://github.com/Bahex/tagore, to be installed. Make sure tagore is present in PATH.

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


Usage and arguments:
```
Rscript FLARE_Visualizer.R --help

Options:
	-i PATH, --input-file=PATH
		path of the ancestry file produced by FLARE. it can either be a vcf file or vcf.gz file.

	-d, --do-parallel
		Perform parallel computing for faster processing. Turned off by default

	-o PATH, --output-dir=PATH
		path of the output directory

	-h, --help
		Show this help message and exit
```