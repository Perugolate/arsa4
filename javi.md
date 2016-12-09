```r
library(readr)
library(magrittr)
library(lubridate)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
df1 <- read_csv("javi.csv")
df1$strain <- df1$line
df1 <- separate(df1, line, into = c("treatment", "line"))
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
```r
unique(df1$plate)
```
```
[1] 1 2 3
```
So he ran three plates.

Which strains did he look at:

```r
unique(df1$line) %>% length
```
``` 
 [1] "WT1C2"       "WT1C3"       "WT2C3"       "WT3C2"       "WT3C3"
 [6] "WT4C2"       "WT4C3"       "WT5C2"       "WT5C3"       "T1 1 C2 L"
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

Javi measured 63 strains (plus "BLANK"). I think this is pretty much all of them.

Did he balance these across plates?:

```r
dplyr::filter(df1, plate == 1) %>% select(line) %>% unique %>% as.data.frame
```
```
        line
1      WT1C2
2      WT1C3
3      WT2C3
4      WT3C2
5      WT3C3
6      WT4C2
7      WT4C3
8      WT5C2
9      WT5C3
10 T1 1 C2 L
11 T1 1 C2 S
12 T1 1 C3 L
13 T1 1 C3 S
14 T1 2 C2 L
15 T1 2 C2 S
16 T1 2 C3 L
17   T1 3 C2
18 T1 3 C2 L
19 T1 3 C2 S
20   T1 3 C3
21     BLANK
22 T1 3 C3 S
```
```r
dplyr::filter(df1, plate == 2) %>% select(line) %>% unique %>% as.data.frame
```
```
        line
1     WT 1 C
2     WT 2 C
3     WT 3 C
4     WT 4 C
5     WT 5 C
6   ANCESTOR
7    WT 2 C2
8   T1 1 L C
9   T1 1 S C
10  T1 2 L C
11  T1 2 S C
12    T1 4 C
13 T1 4 C2 L
14 T1 4 C3 L
15 T1 4 C3 S
16    T1 5 C
17   T1 5 C2
18   T1 5 C3
19  T1T2 1 C
20 T1T2 1 C2
21 T1T2 1 C3
22     BLANK
```
```r
dplyr::filter(df1, plate == 3) %>% select(line) %>% unique %>% as.data.frame
```
```
          line
1     T1T2 2 C
2    T1T2 2 C2
3    T1T2 2 C3
4  T1T2 4 C2 S
5   T1T2 4 C3L
6   T1T2 4 L C
7   T1T2 4 S C
8  T1T2 5 C2 L
9  T1T2 5 C2 S
10 T1T2 5 C3 S
11  T1T2 5 L C
12    T1 3 L C
13    T1 3 S C
14    T1T2 3 C
15   T1T2 3 C2
16   T1T2 3 C3
17 T1T2 4 C2 L
18 T1T2 4 C3 S
19 T1T2 5 C3 L
20  T1T2 5 S C
21   T1 3 C3 L
22       BLANK
```

Did he fuck. He even ran a whole plate with no unselected controls on it. Also, I wonder if these are pseudo replicates ...

Parameters by treatment within line


```r
tre_line <- filter(df1, strain != "BLANK" & strain != "ANCESTOR") %>% group_by(treatment, line) %>% summarise(vmax = mean(vmax), lag = mean(lag), final_OD = mean(final_OD))
```
```r
by_strain <- filter(df1, strain != "BLANK") %>% group_by(strain) %>% summarise(vmax = mean(vmax), lag = mean(lag), final_OD = mean(final_OD)) %>% separate(strain, into = c("treatment", "line"), remove = FALSE)
lag <- ggplot(by_strain, aes(x=treatment, y=lag, color=line)) + geom_point(size=5)
vmax <- ggplot(by_strain, aes(x=treatment, y=vmax, color=line)) + geom_point(size=5)
final_OD <- ggplot(by_strain, aes(x=treatment, y=final_OD, color=line)) + geom_point(size=5)
plot_grid(lag, vmax, final_OD, nrow = 1)
```

