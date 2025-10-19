# Expected utility and global minimum variance shrinkage portfolios.

<!-- badges: start -->
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/HDShOP)](https://cran.r-project.org/package=HDShOP)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/Otryakhin-Dmitry/global-minimum-variance-portfolio/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Otryakhin-Dmitry/global-minimum-variance-portfolio/actions/workflows/R-CMD-check.yaml)
[![](https://cranlogs.r-pkg.org/badges/grand-total/HDShOP?color=orange)](https://cranlogs.r-pkg.org/)
[![](https://cranlogs.r-pkg.org/badges/HDShOP)](https://cranlogs.r-pkg.org/)
<!-- badges: end -->

## Description
The package features a framework for working with high-dimensional shrinkage 
optimal portfolios. It allows constructing those in two ways: 1) by applying 
shrinkage directly to the portfolio weights (function `MVShrinkPortfolio`) and 
2) by obtaining shrinkage estimates of mean returns and covariance matrices 
(function `MeanVar_portfolio`).

## Installation

The latest stable release is always on CRAN:
``` r
install.packages('HDShOP')
```
The latest development version can be installed in the following way:
``` r
library("remotes")

u<-"Otryakhin-Dmitry/"
r<-"global-minimum-variance-portfolio"

re <- paste(u,r,sep="")
remotes::install_github(repo=re, subdir="")
```

## Example
In this example, returns of assets from S&P500 are loaded and an MV portfolio is
created, for which methods `summary` and `plot` are called.  
```r
library(HDShOP)

# loading S&P daily asset returns
data("SP_daily_asset_returns")
assets <- t(SP_daily_asset_returns[2:301, 2:201])

gamma<-1
p <- nrow(assets)
b<-exp(-0.1*(1:p))

# creating an MV shrinkage portfolio
sh_mv_port <- MVShrinkPortfolio(x=assets, gamma=gamma,
                                type='shrinkage', b=b, beta = 0.05)

# Making a summary and plotting the portfolio
summary(sh_mv_port)
plot(sh_mv_port)

```
## Citation
Taras Bodnar, Solomiia Dmytriv, Yarema Okhrin, Dmitry Otryakhin & Nestor
Parolya (12 May 2025): High-Dimensional portfolio selection with HDShOP package, The
European Journal of Finance, [DOI: 10.1080/1351847X.2025.2501637](https://doi.org/10.1080/1351847X.2025.2501637)
