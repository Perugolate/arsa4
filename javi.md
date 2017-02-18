# Table of contents

- [Prep data](#prep-data)
 * [Rescue the IDs](#rescue-the-ids)
 * [Check the design](#check-the-design)
 * [Growth parameters by treatment/line/strain](#growth-parameters-by-treatmentlinestrain)
- [Effect of mutation on growth](#effect-of-mutation-on-growth)
 * [Growth parameters by treatment/mutation](#growth-parameters-by-treatmentmutation)
 * [Models of growth parameters by mutation](#models-of-growth-parameters-by-mutation)
- [Model summaries](#model-summaries)
 * [`vmax ~ mutation`](#vmax--mutation) 
 * [`lag ~ mutation`](#lag--mutation) 
 * [`final_OD ~ mutation`](#final_od--mutation) 

# Prep data

```r
library(readr)
library(magrittr)
library(lubridate)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(visreg)
df1 <- read_csv("javi.csv")
df1$strain <- df1$line
df1 <- separate(df1, line, into = c("treatment", "line"))
# convert lag time into decimal minutes
df1$lag <- hms(df1$lag)
df1$lag <- hour(df1$lag) * 60 + minute(df1$lag) + second(df1$lag) / 60
df1
```

```
# A tibble: 1,152 × 10
   treatment  line   rep  well  vmax     t_vmax      lag final_OD plate  strain
*      <chr> <chr> <int> <chr> <dbl>     <time>    <dbl>    <dbl> <int>   <chr>
1         WT     1     1    A1   4.0  8700 secs 50.75000    1.096     1 WT 1 C2
2         WT     1     2    B1   4.0  9300 secs 54.25000    1.070     1 WT 1 C2
3         WT     1     3    C1   3.9  9300 secs 51.53333    1.071     1 WT 1 C2
4         WT     1     4    D1   3.9  8700 secs 52.56667    1.078     1 WT 1 C2
5         WT     1     5    E1   4.0  9900 secs 55.75000    1.073     1 WT 1 C2
6         WT     1     6    F1   4.0  9300 secs 56.75000    1.087     1 WT 1 C2
7         WT     1     7    G1   4.0 10500 secs 56.50000    1.080     1 WT 1 C2
8         WT     1     8    H1   3.9  9900 secs 54.36667    1.095     1 WT 1 C2
9         WT     1     9    I1   4.0  9300 secs 54.50000    1.087     1 WT 1 C2
10        WT     1    10    J1   4.1  9300 secs 57.31667    1.096     1 WT 1 C2
# ... with 1,142 more rows
```

## Rescue the IDs

The strain names are formatted differently from the previous data. Rescue them with the unique ID numbers.

```r
ids <- read_csv("ids_javi_to_olga.csv")
df2 <- left_join(df1, ids, by = c("strain" = "javi_id"))
```

## Check the design

```r
unique(df2$plate)
```

```
[1] 1 2 3
```

So there were three plates.

Which strains were on the plates?:

```r
unique(df2$strain)
```

``` 
 [1] "WT 1 C2"     "WT 1 C3"     "WT 2 C3"     "WT 3 C2"     "WT 3 C3"
 [6] "WT 4 C2"     "WT 4 C3"     "WT 5 C2"     "WT 5 C3"     "T1 1 C2 L"
[11] "T1 1 C2 S"   "T1 1 C3 L"   "T1 1 C3 S"   "T1 2 C2 L"   "T1 2 C2 S"
[16] "T1 2 C3 L"   "T1 3 C2"     "T1 3 C2 L"   "T1 3 C2 S"   "T1 3 C3"
[21] "BLANK"       "T1 3 C3 S"   "WT 1 C"      "WT 2 C"      "WT 3 C"
[26] "WT 4 C"      "WT 5 C"      "ANCESTOR"    "WT 2 C2"     "T1 1 L C"
[31] "T1 1 S C"    "T1 2 L C"    "T1 2 S C"    "T1 4 C"      "T1 4 C2 L"
[36] "T1 4 C3 L"   "T1 4 C3 S"   "T1 5 C"      "T1 5 C2"     "T1 5 C3"
[41] "T1T2 1 C"    "T1T2 1 C2"   "T1T2 1 C3"   "T1T2 2 C"    "T1T2 2 C2"
[46] "T1T2 2 C3"   "T1T2 4 C2 S" "T1T2 4 C3L"  "T1T2 4 L C"  "T1T2 4 S C"
[51] "T1T2 5 C2 L" "T1T2 5 C2 S" "T1T2 5 C3 S" "T1T2 5 L C"  "T1 3 L C"
[56] "T1 3 S C"    "T1T2 3 C"    "T1T2 3 C2"   "T1T2 3 C3"   "T1T2 4 C2 L"
[61] "T1T2 4 C3 S" "T1T2 5 C3 L" "T1T2 5 S C"  "T1 3 C3 L"
```

```r
unique(df2$strain) %>% length
```

```
[1] 64
```

So 63 strains were measured (plus "BLANK"). There are three strains missing (T1-2-C3-S, T1-4-C2-S, and T1-4-C2-S).

Is everything balanced across plates?:

```r
dplyr::filter(df2, plate == 1) %>% select(strain) %>% table
```

```
    BLANK T1 1 C2 L T1 1 C2 S T1 1 C3 L T1 1 C3 S T1 2 C2 L T1 2 C2 S T1 2 C3 L
       48        16        16        16        16        16        16        16
  T1 3 C2 T1 3 C2 L T1 3 C2 S   T1 3 C3 T1 3 C3 S   WT 1 C2   WT 1 C3   WT 2 C3
       16        16        16        16        16        16        16        16
  WT 3 C2   WT 3 C3   WT 4 C2   WT 4 C3   WT 5 C2   WT 5 C3
       16        16        16        16        16        16
```

```r
dplyr::filter(df2, plate == 2) %>% select(strain) %>% table
```

```
    BLANK T1 1 C2 L T1 1 C2 S T1 1 C3 L T1 1 C3 S T1 2 C2 L T1 2 C2 S T1 2 C3 L
       48        16        16        16        16        16        16        16
  T1 3 C2 T1 3 C2 L T1 3 C2 S   T1 3 C3 T1 3 C3 S   WT 1 C2   WT 1 C3   WT 2 C3
       16        16        16        16        16        16        16        16
  WT 3 C2   WT 3 C3   WT 4 C2   WT 4 C3   WT 5 C2   WT 5 C3
       16        16        16        16        16        16
```

```r
dplyr::filter(df2, plate == 3) %>% select(strain) %>% table
```

```
      BLANK   T1 3 C3 L    T1 3 L C    T1 3 S C    T1T2 2 C   T1T2 2 C2
         48          16          16          16          16          16
  T1T2 2 C3    T1T2 3 C   T1T2 3 C2   T1T2 3 C3 T1T2 4 C2 L T1T2 4 C2 S
         16          16          16          16          16          16
 T1T2 4 C3L T1T2 4 C3 S  T1T2 4 L C  T1T2 4 S C T1T2 5 C2 L T1T2 5 C2 S
         16          16          16          16          16          16
T1T2 5 C3 L T1T2 5 C3 S  T1T2 5 L C  T1T2 5 S C
         16          16          16          16
```

Nope. Plate three does not have any unselected controls on it. Treatment T1T2 is also only present on plate three. Also, I wonder if these 16 replicates are pseudo replicates ...

## Growth parameters by treatment/line/strain

```r
by_strain <- filter(df2, strain != "BLANK") %>% group_by(strain) %>%
  summarise(vmax = mean(vmax), lag = mean(lag), final_OD = mean(final_OD)) %>%
  separate(strain, into = c("treatment", "line"), remove = FALSE)
lag <- ggplot(by_strain, aes(x = treatment, y = lag, color = line)) +
  geom_point(size = 5)
vmax <- ggplot(by_strain, aes(x = treatment, y = vmax, color = line)) +
  geom_point(size = 5)
final_OD <- ggplot(by_strain, aes(x = treatment, y = final_OD, color = line)) +
  geom_point(size = 5)
plot_grid(lag, vmax, final_OD, nrow = 1)
# summarise by ID in order to join with mutation data
by_id <- filter(df2, strain != "BLANK") %>% group_by(ID) %>%
  summarise(vmax = mean(vmax), lag = mean(lag), final_OD = mean(final_OD))
```

![](https://github.com/Perugolate/arsa4/blob/master/plots/growth_by_strains.png)

![](https://github.com/Perugolate/arsa4/blob/master/plots/growth_by_tre_box.png)

# Effect of mutation on growth

These are just quick and dirty models incorporating the two most frequently mutated operons (ytr and gra - all resistant strains except one have a mutation in one of these two operons). i.e. agglomerating at the operon level - I will fit a proper model at some point. Perhaps something like `vmax ~ ytrA + ytrB + graR + graS + vraS + tcaA ...` where each variable is a factor and the levels are the genotypes.

Combine with previous mutation data and test for effect of ytr/gra mutations on growth rate:

```r
odf1 <- read.csv("o90.csv", skip = 11) # skip the comment lines
# line combines multiple factors so separate to get treatment on its own.
odf2 <- separate(odf1, line, into=c("treatment", "line", "colony"), sep = "-",
  extra = "merge")
# get rid of the populations since we don't have growth measurements for them
odf4 <- subset(odf2, colony != "P") %>% droplevels
# agglomerate mutations for ytrAB, and graRS
odf7 <- data.frame(dplyr::select(odf4, ID, treatment, line),
  ytr = rowSums(select(odf4, ytrA, ytrB)), gra = rowSums(select(odf4, graS, graR)))
# join with Javi's growth data
gro_mu <- left_join(by_id, odf7)
## this part is super hacky
# gather into long
gro_mu_l <- gather(gro_mu, key = mutation, value = presence, ytr, gra)
# keep the strains which have a mutation
foo <- dplyr::filter(gro_mu_l, presence == 1)
# create an object with only the controls
bar <- filter(gro_mu_l, treatment == "con" & mutation == "ytr")
# set ehe level to none
bar$mutation <- gsub("ytr", "none", bar$mutation)
# combine to make an object with variable "mutation" with values "none", "ytr", and "gra".
groMU <- rbind(foo, bar)
```

## Growth parameters by treatment/mutation

```r
png("plots/growth_by_tremu.png", height = 480 * 0.8, width = 3*(480 * 0.8))
mvmax <- ggplot(groMU, aes(x = treatment, y = vmax, color = mutation)) +
  geom_point(size = 5, position = position_dodge(width = 0.5))
mlag <- ggplot(groMU, aes(x = treatment, y = lag, color = mutation)) +
  geom_point(size = 5, position=position_dodge(width = 0.5))
mod <- ggplot(groMU, aes(x = treatment, y = final_OD, color = mutation)) +
  geom_point(size = 5, position=position_dodge(width = 0.5))
plot_grid(mvmax, mlag, mod, nrow = 1)
dev.off()
```

![](https://github.com/Perugolate/arsa4/blob/master/plots/growth_by_tremu.png)

## Models of growth parameters by mutation

```r
png("plots/growth_by_mu.png", height = 480 * 0.8, width = 3 * (480 * 0.8))
par(mfrow = c(1,3), cex = 1.2)
glm(vmax ~ mutation, data=groMU, family = gaussian) %>%
  visreg(main = "vmax", ylab = "vmax")
glm(lag ~ mutation, data=groMU, family = gaussian) %>%
  visreg(main = "lag", ylab = "lag phase (minutes)")
glm(final_OD ~ mutation, data=groMU, family = gaussian) %>%
  visreg(main = "final OD", ylab = "final OD (600 nm)")
dev.off()
```

![](https://github.com/Perugolate/arsa4/blob/master/plots/growth_by_mu.png)


![](https://github.com/Perugolate/arsa4/blob/master/plots/grwoth_by_mu_plus_tre.png)

# Model summaries

## `vmax ~ mutation`

```r
glm(vmax ~ mutation, data = groMU, family = gaussian) %>% summary
```

```
Call:
glm(formula = vmax ~ mutation, family = gaussian, data = groMU)

Deviance Residuals:
     Min        1Q    Median        3Q       Max
-1.09407  -0.16282   0.02333   0.36843   0.75593

Coefficients:
             Estimate Std. Error t value Pr(>|t|)
(Intercept)    2.7114     0.1026  26.430  < 2e-16 ***
mutationnone   1.2340     0.1498   8.236 2.47e-11 ***
mutationytr    0.3639     0.1292   2.817  0.00662 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for gaussian family taken to be 0.1789134)

    Null deviance: 23.212  on 60  degrees of freedom
Residual deviance: 10.377  on 58  degrees of freedom
AIC: 73.062

Number of Fisher Scoring iterations: 2
```

## `lag ~ mutation`

```r
glm(lag ~ mutation, data = groMU, family = gaussian) %>% summary
```

```
Call:
glm(formula = lag ~ mutation, family = gaussian, data = groMU)

Deviance Residuals:
    Min       1Q   Median       3Q      Max
-59.098  -21.619   -1.197   15.464   81.967

Coefficients:
             Estimate Std. Error t value Pr(>|t|)
(Intercept)   135.966      6.979  19.481  < 2e-16 ***
mutationnone  -69.539     10.194  -6.822 5.84e-09 ***
mutationytr     1.202      8.790   0.137    0.892
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for gaussian family taken to be 828.0824)

    Null deviance: 103942  on 60  degrees of freedom
Residual deviance:  48029  on 58  degrees of freedom
AIC: 587.9

Number of Fisher Scoring iterations: 2
```

## `final_OD ~ mutation`

```r
glm(final_OD ~ mutation, data = groMU, family = gaussian) %>% summary
```

```
Call:
glm(formula = final_OD ~ mutation, family = gaussian, data = groMU)

Deviance Residuals:
      Min         1Q     Median         3Q        Max
-0.202439  -0.056064   0.000904   0.050811   0.202923

Coefficients:
             Estimate Std. Error t value Pr(>|t|)
(Intercept)   0.94883    0.01963  48.338  < 2e-16 ***
mutationnone  0.09027    0.02867   3.149 0.002593 **
mutationytr   0.10049    0.02472   4.065 0.000147 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for gaussian family taken to be 0.006550187)

    Null deviance: 0.49633  on 60  degrees of freedom
Residual deviance: 0.37991  on 58  degrees of freedom
AIC: -128.69

Number of Fisher Scoring iterations: 2
```

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
