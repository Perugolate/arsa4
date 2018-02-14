## Proper model

Need to alter this to work with merged.vcf (test.vcf is derived from `snippy-core` which excludes indels).

Could also check output from `cat merged.vcf | vcf-to-tab > merged.tab`.

```r
library(vcfR) # install via devtools::install_github(repo="knausb/vcfR")
library(magrittr)
library(tidyr)
library(dplyr)
tf1 <- read.vcfR("test.vcf") %>% vcfR2tidy
tf1 <- tf1$gt
tf2 <- separate(tf1, Indiv, into = "ID", extra = "drop") %>% filter(ID != "Reference")
# there are no SNPs with more than 1 alternative allele
# try joining this to summarized growth data - remember ID is chr here
select(tf2, POS, ID, gt_GT) %>% group_by(POS, ID) %>% summarize(gt_GT) %>% spread(POS, gt_GT)
```

Some visualizations

```r
# Find the files.
vcf_file <- "merged.vcf" # derived this from vcf-merge on *.filt.subs.vcf.gz from snippy
dna_file <- "sh1000_pol_6.fa"
gff_file <- "sh1000_pol_6.gff"

# Input the files.
vcf <- read.vcfR(vcf_file, verbose = FALSE)
dna <- ape::read.dna(dna_file, format = "fasta")
gff <- read.table(gff_file, sep="\t", quote="")

# Create a chromR object.
chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=dna, ann=gff, verbose=TRUE)
#chrom <- masker(chrom, min_QUAL=0, min_DP=350, max_DP=650, min_MQ=59.5, max_MQ=60.5)
chrom <- proc.chromR(chrom, verbose = TRUE)
#chrom <- proc.chromR(chrom, verbose=FALSE, win.size=1e3)
chromoqc(chrom)
#plot(chrom)
extract.gt(chrom, as.numeric = TRUE) %>% heatmap.bp
```

Might try `snpeff` on the merge vcf:

```sh
mkdir filt && cd filt
# get all the filtered vcfs from snippy output
cp ../[0-9]*_S*[0-9]/*_S*filt.subs.vcf.gz* .
vcf-merge *vcf.gz > merged.vcf
# add an ANN field describing the predicted snpeff
~/opt/snippy-3.1/binaries/noarch/snpEff ann -no-downstream -no-upstream -no-intergenic -no-utr -c reference/snpeff.config -dataDir . -noStats ref merged.vcf > eff.vcf
# convert to tab format (might be easier to wrangle for fitting a model)
~/opt/snippy-3.1/bin/snippy-vcf_to_tab --gff reference/ref.gff --ref reference/ref.fa --vcf eff.vcf > eff.tab
```

OK, so the result of the above chunk is that the sample info is stripped from the output of `snippy-vcf_to_tab`. snpeff handles the multi vcf OK so it contains ANN. Need a way of converting it to tab or I just go witht the vcf. Look at code for `snippy-vcf_to_tab` since this drops the columns. Might also be possible to merge the individual tab files.

with `eff.vcf` and `eff.tab`, try using the vcf then joing to the tab to rescue ANN:
```r
eff1 <- read.vcfR("eff.vcf")
eff2 <- vcfR2tidy(eff1)
eff3 <- eff2$gt
eff4 <- select(eff3, POS, Indiv, gt_GT_alleles)
# read in the tabix file to rescue the snpeff annotations
tab1 <- read_tsv("eff.tab")
# try combining gene with pos since the intergenic mutations are NA
tab2 <- unite(tab1, locus, c(GENE, POS), remove = FALSE)
df7 <- left_join(eff4, tab2)
df8 <- select(df7, Indiv, locus, gt_GT_alleles)
df9 <- spread(df8, locus, gt_GT_alleles)
df9[is.na(df9)] <- "Z"
df9 <- separate(df9, Indiv, into = c("ID", "sample"))
# join with by_id and fit a model
by_id$ID <- as.character(by_id$ID)
foo <- left_join(by_id, df9)
```

```r
df3 <- select(df2, treatment, ID, line)
df3$ID <- as.character(df3$ID)
df4 <- left_join(foo, unique(df3))
df4 <- select(df4, line, ID, treatment, everything())
ggplot(df4, aes(x = treatment, y = vmax, color = ytrA_1935783)) +
  geom_point(size = 5, position = position_dodge(width = 0.5))
```
`
