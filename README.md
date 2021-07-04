README
================
Michael B. Sohn
7/01/2021

## CMMB: Compositional Mediation Model for Binary Outcomes

The cmmb function estimates direct and indirect effects of treatment on
binary outcomes mediated through a compositional mediator that consists
of multiple components. For detailed information about the arguments,
please see the documentation for *cmmb()*.

Note: the number of components can be high-dimensional. However, it will
be computationally intensive to run it with high-dimensional components
as it tests the indirect effect using debiased bootstrap estimates based
on 2,000 (in default) random samplings.

### Installation

install.packages(“devtools”)

devtools::install\_github(“mbsohn/cmmb”)

### Example: estimate direct and total indirect effects

``` r
### load functions used in data generation and performance comparison
library(cmmb)
source("comparators.R")
### Generate a simulated dataset of 50 samples and 5 taxa
set.seed(2021)
### The gen.cmm.sim.data function simulates data based on the logistic normal distribution. The first argument takes the number of samples and the second takes the number of components (taxa).
sim.dat <- gen.cmm.sim.data(50, 5)
### Run CMM for binary outcomes
rslt <- cmmb(Y=sim.dat$Y, M=sim.dat$M, tr=sim.dat$tr, X=sim.dat$X)
rslt
```

    ## $total
    ##       Estimate Lower Limit Upper Limit
    ## DE   0.1831665   0.1229387   0.2930580
    ## TIDE 0.2585639   0.1391205   0.3345654
    ## 
    ## $cwprod
    ##           Estimate Lower Limit Upper Limit
    ## taxon1  0.45250088  0.06145777   0.7646973
    ## taxon2  0.00179244 -0.18731449   0.1584444
    ## taxon3 -0.20841877 -0.53383632   0.1692616
    ## taxon4  0.61865737 -0.14664985   1.1418571
    ## taxon5  0.51405582 -0.45317917   1.1773031
    ## 
    ## attr(,"class")
    ## [1] "cmmb"

``` r
### Plot products of component-wise path coefficients
plot_cw_ide(rslt)
```

![](README_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

### Example: sensitivity analysis for the total indirect effect

``` r
rslt.sa <- cmmb(Y=sim.dat$Y, M=sim.dat$M, tr=sim.dat$tr, X=sim.dat$X, ForSA=TRUE)
### Plot sensitivity of the estimated total indirect effect
plot_cmmb_sa(rslt.sa)
```

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

### Performance comparison

``` r
set.seed(2021)
sim.dat <- gen.cmm.sim.data(50, 5)
PCS(Y=sim.dat$Y, M=sim.dat$M, tr=sim.dat$tr, X=sim.dat$X)
```

    ##           Estimate Lower Limit Upper Limit
    ## DE   -0.0006223289   -0.573698   0.8385662
    ## TIDE  0.1600846113   -2.155939   2.4267162

``` r
PCP(Y=sim.dat$Y, M=sim.dat$M, tr=sim.dat$tr, X=sim.dat$X)
```

    ##       Estimate Lower Limit Upper Limit
    ## DE   0.4499679   0.2390003   0.7423284
    ## TIDE 0.0217547  -0.2973603   0.1997812

``` r
cmmb(Y=sim.dat$Y, M=sim.dat$M, tr=sim.dat$tr, X=sim.dat$X)$total
```

    ##       Estimate Lower Limit Upper Limit
    ## DE   0.1831665   0.1216986   0.2847858
    ## TIDE 0.2585639   0.1530183   0.3294128
