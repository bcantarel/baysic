# baysic

The basics script requires:
[vcftools](https://vcftools.github.io/index.html)
R packages:
[R2jags](https://cran.r-project.org/web/packages/R2jags/index.html)
[getopt](https://cran.r-project.org/web/packages/getopt/index.html)
[reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)

With these things installed and available in your path, you can run BAYSIC using:

baysic_blc.pl -f reference.fasta -p output_prefix -l /full/path/lca.R *.vcf  

where the reference.fa is the reference for which the reads were aligned and output_prefix will name your output files
