- [Dependencies](#dependencies)
- [Mutation data](#mutation-data)
  - [Number of mutations by treatment](#number-of-mutations-by-treatment)
- [Growth data](#growth-data)
  - [Growth parameters by treatment](#growth-parameters-by-treatment)
  - [Growth parameters by treatment and mutation](#growth-parameters-by-treatment-and-mutation)
- [Models of growth parameters by mutation](#models-of-growth-parameters-by-mutation)
  - [Model summaries](#model-summaries)

# Dependencies
```r
library(readr)
library(tidyr)
library(magrittr)
library(visreg)
library(contrast)
library(dplyr)
library(plotrix)
library(ggplot2)
library(cowplot)
library(nlme)
library(lubridate)
```

# Mutation data

```r
odf1 <- read.csv("o90.csv", skip = 11) # skip the comment lines
# line combines multiple factors so separate to get treatment on its own.
odf2 <- separate(odf1, line, into = c("treatment", "line", "colony"),
		sep = "-", extra = "merge")
# exclude populations - less power to detect variants therefore can't be in the same model as colonies.
odf4 <- subset(odf2, colony != "P") %>% droplevels
odf5 <- data.frame(odf4$treatment, odf4$line,
		   rowSums(dplyr::select(odf4, -ID, -treatment, -line, -colony))) %>% droplevels
colnames(odf5) <- c("treatment", "line", "mutations")
# summarise by treatment for plot
odf6 <- group_by(odf5, treatment)
odf6p <- summarise(odf6, m_mu = mean(mutations), ci = 1.96 * std.error(mutations))
# contrast T1 with T1T2, line as random factor
mixm <- lme(mutations ~ treatment, random = ~1|line, data = odf5, method = "REML")
contrast(mixm, list(treatment = "T1"), list(treatment = "T1T2"))
```
```
lme model parameter contrast

    Contrast      S.E.      Lower     Upper     t df Pr(>|t|)
1 -0.1945512 0.2048562 -0.6043247 0.2152222 -0.95 60   0.3461
```
```r
summary(mixm)
```
```
Linear mixed-effects model fit by REML
 Data: odf5
       AIC      BIC    logLik
  160.6369 171.2725 -75.31843

Random effects:
 Formula: ~1 | line
        (Intercept)  Residual
StdDev:   0.7352224 0.6937649

Fixed effects: mutations ~ treatment
                 Value Std.Error DF  t-value p-value
(Intercept)   1.533333 0.3744299 58 4.095115   1e-04
treatmentT1   1.757946 0.2219247 58 7.921364   0e+00
treatmentT1T2 1.952498 0.2364301 58 8.258245   0e+00
 Correlation:
              (Intr) trtmT1
treatmentT1   -0.386
treatmentT1T2 -0.362  0.602

Standardized Within-Group Residuals:
       Min         Q1        Med         Q3        Max
-2.5981418 -0.5774465  0.1724957  0.5835358  2.2297958

Number of Observations: 65
Number of Groups: 5
```

## Number of mutations by treatment

```r
png("plots/n_mu_by_tre.png", width = 2*480)
gpoint <- ggplot(odf6p, aes(x = treatment, y = m_mu)) + geom_point(size = 5) +
  geom_errorbar(aes(ymin = m_mu - ci, ymax = m_mu + ci), width = 0.2) +
  ylab("Number of mutations")
gbox <- ggplot(odf5, aes(x = treatment, y = mutations)) + geom_boxplot() +
  ylab("Number of mutations")
plot_grid(gpoint, gbox)
dev.off()
```

![](https://github.com/Perugolate/arsa4/blob/master/plots/n_mu_by_tre.png)

```r
png("plots/n_mu_by_tre_line.png")
lm(mutations ~ treatment/line, data = odf5) %>% visreg("treatment", by = "line")
dev.off()
```

![](https://github.com/Perugolate/arsa4/blob/master/plots/n_mu_by_tre_line.png)

# Growth data

```r
df1 <- read_csv("javi.csv")
df1$strain <- df1$line
df1 <- separate(df1, line, into = c("treatment", "line"))
# convert lag time into decimal minutes
df1$lag <- hms(df1$lag)
df1$lag <- hour(df1$lag) * 60 + minute(df1$lag) + second(df1$lag) / 60
# rescue strain names with the unique ID numbers
ids <- read_csv("ids_javi_to_olga.csv")
df2 <- left_join(df1, ids, by = c("strain" = "javi_id"))
# summarise by strain
by_strain <- filter(df2, strain != "BLANK") %>% group_by(strain) %>%
  summarise(vmax = mean(vmax), lag = mean(lag), final_OD = mean(final_OD)) %>%
  separate(strain, into = c("treatment", "line"), remove = FALSE)
# summarise by ID in order to join with mutation data
by_id <- filter(df2, strain != "BLANK") %>% group_by(ID) %>%
  summarise(vmax = mean(vmax), lag = mean(lag), final_OD = mean(final_OD), plate = unique(plate))
# combine with previous mutation data
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
# set the level to none
bar$mutation <- gsub("ytr", "none", bar$mutation)
# combine to make an object with variable "mutation" with values "none", "ytr", and "gra".
groMU <- rbind(foo, bar)
groMU$mutation <- as.factor(groMU$mutation)
groMU$mutation <- relevel(groMU$mutation, ref = "none")
```

## Growth parameters by treatment

```r
png("plots/growth_by_tre_box2.png", width = 3*480)
mvmax <- ggplot(groMU, aes(x = treatment, y = vmax)) + geom_boxplot() +
  ylab(expression(V["max"]))
mlag <- ggplot(groMU, aes(x = treatment, y = lag)) + geom_boxplot() +
  ylab("lag phase (minutes)")
mod <- ggplot(groMU, aes(x = treatment, y = final_OD)) + geom_boxplot() +
  ylab(expression(paste("final ", "OD"["600"], sep = " ")))
plot_grid(mvmax, mlag, mod, nrow = 1)
dev.off()
```

![](https://github.com/Perugolate/arsa4/blob/master/plots/growth_by_tre_box2.png)

## Growth parameters by treatment and mutation

```r
png("plots/growth_by_tremu.png", height = 480 * 0.8, width = 3*(480 * 0.8))
mvmax <- ggplot(groMU, aes(x = treatment, y = vmax, color = mutation)) +
  geom_point(size = 5, position = position_dodge(width = 0.5)) +
  ylab(expression(V["max"]))
mlag <- ggplot(groMU, aes(x = treatment, y = lag, color = mutation)) +
  geom_point(size = 5, position=position_dodge(width = 0.5)) +
  ylab("lag phase (minutes)")
mod <- ggplot(groMU, aes(x = treatment, y = final_OD, color = mutation)) +
  geom_point(size = 5, position=position_dodge(width = 0.5)) +
  ylab(expression(paste("final ", "OD"["600"], sep = " ")))
plot_grid(mvmax, mlag, mod, nrow = 1)
dev.off()
```

![](https://github.com/Perugolate/arsa4/blob/master/plots/growth_by_tremu.png)

# Models of growth parameters by mutation

```r
png("plots/growth_by_mu.png", height = 480 * 0.8, width = 3 * (480 * 0.8))
par(mfrow = c(1,3), cex = 1.2)
glm(vmax ~ mutation, data=groMU, family = gaussian) %>%
  visreg(main = expression(V["max"]), ylab = expression(V["max"]))
glm(lag ~ mutation, data=groMU, family = gaussian) %>%
  visreg(main = "lag phase", ylab = "lag phase (minutes)")
glm(final_OD ~ mutation, data=groMU, family = gaussian) %>%
  visreg(main = "final OD", ylab = expression(paste("final ", "OD"["600"])))
dev.off()
```

![](https://github.com/Perugolate/arsa4/blob/master/plots/growth_by_mu.png)

```r
png("plots/growth_by_mu_plus_tre.png", height = 480 * 0.8, width = 3*(480 * 0.8))
gvmax <- ggplot(groMU, aes(x = mutation, y = vmax, color = treatment)) + 
  geom_point(size = 5, position=position_dodge(width = 0.5)) +
  ylab(expression(V["max"]))
glag <- ggplot(groMU, aes(x = mutation, y = lag, color = treatment)) + 
  geom_point(size = 5, position=position_dodge(width = 0.5)) +
  ylab("lag phase (minutes)")
glod <- ggplot(groMU, aes(x = mutation, y = final_OD, color = treatment)) + 
  geom_point(size = 5, position=position_dodge(width = 0.5)) +
  ylab(expression(paste("final ", "OD"["600"], sep = " ")))
plot_grid(gvmax, glag, glod, nrow = 1)
dev.off()
```

![](https://github.com/Perugolate/arsa4/blob/master/plots/growth_by_mu_plus_tre.png)

## Model summaries

### `vmax ~ mutation`

```r
lme(vmax ~ mutation, random = ~1|line, data = groMU, method = "REML") %>% summary
```

```
Linear mixed-effects model fit by REML
 Data: groMU
       AIC     BIC    logLik
  82.23238 92.5346 -36.11619

Random effects:
 Formula: ~1 | line
        (Intercept)  Residual
StdDev:   0.1338896 0.4070154

Fixed effects: vmax ~ mutation
                Value Std.Error DF  t-value p-value
(Intercept)  3.945417 0.1209520 54 32.61969       0
mutationgra -1.255830 0.1539462 54 -8.15759       0
mutationytr -0.854967 0.1324057 54 -6.45717       0
 Correlation:
            (Intr) mttngr
mutationgra -0.593
mutationytr -0.690  0.470

Standardized Within-Group Residuals:
       Min         Q1        Med         Q3        Max
-2.5351730 -0.5588148  0.1481167  0.7649218  1.5152300

Number of Observations: 61
Number of Groups: 5
```

```r
vmax_lme <- lme(vmax ~ mutation, random = ~1|line, data = groMU, method = "REML")
contrast(vmax_lme, list(mutation = "gra"), list(mutation = "ytr"))
```

```
lme model parameter contrast

   Contrast      S.E.     Lower     Upper    t df Pr(>|t|)
1 -0.400863 0.1485144 -0.698373 -0.103353 -2.7 56   0.0092
```

### `lag ~ mutation`

```r
lme(lag ~ mutation, random = ~1|line, data = groMU, method = "REML") %>% summary
```
```
Linear mixed-effects model fit by REML
 Data: groMU
       AIC     BIC    logLik
  550.1268 560.429 -270.0634

Random effects:
 Formula: ~1 | line
        (Intercept) Residual
StdDev:    21.05472 21.74158

Fixed effects: lag ~ mutation
               Value Std.Error DF  t-value p-value
(Intercept) 66.42632 10.962359 54 6.059491       0
mutationgra 62.82264  8.738307 54 7.189338       0
mutationytr 71.26607  7.235584 54 9.849387       0
 Correlation:
            (Intr) mttngr
mutationgra -0.329
mutationytr -0.397  0.363

Standardized Within-Group Residuals:
        Min          Q1         Med          Q3         Max
-2.02508833 -0.54667163 -0.03432731  0.58440504  2.83819265

Number of Observations: 61
Number of Groups: 5
```

```r
lag_lme <- lme(lag ~ mutation, random = ~1|line, data = groMU, method = "REML")
contrast(lag_lme, list(mutation = "gra"), list(mutation = "ytr"))
```
```
lme model parameter contrast

   Contrast     S.E.     Lower   Upper     t df Pr(>|t|)
1 -8.443423 9.102907 -26.67874 9.79189 -0.93 56   0.3576
```

### `final_OD ~ mutation`

```r
lme(final_OD ~ mutation, random = ~1|line, data = groMU, method = "REML") %>% summary
```
```
Linear mixed-effects model fit by REML
 Data: groMU
        AIC       BIC   logLik
  -114.8988 -104.5966 62.44938

Random effects:
 Formula: ~1 | line
        (Intercept)   Residual
StdDev:  0.03919988 0.07280155

Fixed effects: final_OD ~ mutation
                 Value  Std.Error DF  t-value p-value
(Intercept)  1.0390958 0.02570338 54 40.42642  0.0000
mutationgra -0.0903432 0.02848323 54 -3.17180  0.0025
mutationytr  0.0061345 0.02397992 54  0.25582  0.7991
 Correlation:
            (Intr) mttngr
mutationgra -0.483
mutationytr -0.573  0.410

Standardized Within-Group Residuals:
        Min          Q1         Med          Q3         Max
-2.45948392 -0.62907326 -0.06289563  0.72735393  2.62935481

Number of Observations: 61
Number of Groups: 5
```
```r
fod_lme <- lme(final_OD ~ mutation, random = ~1|line, data = groMU, method = "REML")
contrast(fod_lme, list(mutation = "gra"), list(mutation = "ytr"))
```
```
lme model parameter contrast

     Contrast       S.E.      Lower       Upper     t df Pr(>|t|)
1 -0.09647773 0.02875101 -0.1540729 -0.03888253 -3.36 56   0.0014
```


