# Expected utility and global minimum variance portfolios.

<!-- badges: start -->
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/HDShOP)](https://cran.r-project.org/package=HDShOP)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/Otryakhin-Dmitry/global-minimum-variance-portfolio/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Otryakhin-Dmitry/global-minimum-variance-portfolio/actions/workflows/R-CMD-check.yaml)
[![](https://cranlogs.r-pkg.org/badges/grand-total/HDShOP?color=orange)](https://cranlogs.r-pkg.org/)
[![](https://cranlogs.r-pkg.org/badges/HDShOP)](https://cranlogs.r-pkg.org/)
<!-- badges: end -->


## Installation

### The latest version on CRAN:
``` r
install.packages('HDShOP')
```
### The latest development version:
``` r
library("remotes")

u<-"Otryakhin-Dmitry/"
r<-"global-minimum-variance-portfolio"
sd<-"gmvp"

re <- paste(u,r,sep="")
remotes::install_github(repo=re, subdir=sd)
```
