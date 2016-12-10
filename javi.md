# Table of contents

- [Prep data](#prep-data)
 * [Rescue the IDs](#rescue-the-ids)
 * [Check the design](#check-the-design)
- [Growth parameters by treatment/line/strain](#growth-parameters-by-treatmentlinestrain)
- [Effect of mutation on growth](#effect-of-mutation-on-growth)

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

So he ran three plates.

Which strains did he look at:

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
unique(df1$line) %>% length
```

```
[1] 64
```

Javi measured 63 strains (plus "BLANK"). There are three strains missing (T1-2-C3-S, T1-4-C2-S, and T1-4-C2-S).

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

Nope. There is even a whole plate with no unselected controls on it. Also, I wonder if these 16 replicates are pseudo replicates ...

# Growth parameters by treatment/line/strain

```r
by_strain <- filter(df2, strain != "BLANK") %>% group_by(strain) %>% summarise(vmax = mean(vmax), lag = mean(lag), final_OD = mean(final_OD)) %>% separate(strain, into = c("treatment", "line"), remove = FALSE)
lag <- ggplot(by_strain, aes(x=treatment, y=lag, color=line)) + geom_point(size=5)
vmax <- ggplot(by_strain, aes(x=treatment, y=vmax, color=line)) + geom_point(size=5)
final_OD <- ggplot(by_strain, aes(x=treatment, y=final_OD, color=line)) + geom_point(size=5)
plot_grid(lag, vmax, final_OD, nrow = 1)
# the following should be fine for joining with the mutation data
by_id <- filter(df2, strain != "BLANK") %>% group_by(ID) %>% summarise(vmax = mean(vmax), lag = mean(lag), final_OD = mean(final_OD))
```

![](https://github.com/Perugolate/arsa4/blob/master/plots/growth_by_strains.png)

# Effect of mutation on growth

This is just a quick and dirty model incorporating the two most frequently mutated operons (ytr and gra). i.e. agglomerating at the operon level - I will fit a proper model at some point. Perhaps something like `vmax ~ ytrA + ytrB + graR + graS + vraS + tcaA ...` where each variable is a factor and the levels are the genotypes.

Combine with previous mutation data and test for effect of ytr/gra mutations on growth rate:

```r
odf1 <- read.csv("o90.csv", skip=11) # skip the comment lines
# line combines multiple factors so separate to get treatment on its own.
odf2 <- separate(odf1, line, into=c("treatment", "line", "colony"), sep="-", extra="merge")
# get rid of the populations since we don't have growth measurements for them
odf4 <- subset(odf2, colony != "P") %>% droplevels
# agglomerate mutations for ytrAB, and graRS
odf7 <- data.frame(dplyr::select(odf4, ID, treatment, line), ytr=rowSums(select(odf4, ytrA, ytrB)), gra=rowSums(select(odf4, graS, graR)))
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

```r
png("growth_by_mu.png", height = 480*0.8, width = 3*(480*0.8))
par(mfrow = c(1,3))
glm(vmax~mutation, data=groMU, family = gaussian) %>% visreg(main = "vmax")
glm(lag~mutation, data=groMU, family = gaussian) %>% visreg(main = "lag")
glm(final_OD~mutation, data=groMU, family = gaussian) %>% visreg(main = "final OD")
dev.off()
```

![](https://github.com/Perugolate/arsa4/blob/master/plots/growth_by_mu.png)

```r
glm(vmax~mutation, data=groMU, family = gaussian) %>% summary
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
