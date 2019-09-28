# LBRAT

LBRAT is an R program that performs retrospective association testing for longitudinal binary traits in population samples. 


## Installation Instructions:

### Required software and packages
    
1. R (http://www.r-project.org/)
    
2. Package nleqslv, lme4, geepack, mvtnorm, data.table, plyr
    
3. PLINK 1.0 or 1.9 (https://www.cog-genomics.org/plink2)

Please install the required R package before you install LBRAT package. Please install the **LBRAT** as following steps.

 
### Install LBRAT on LUNIX or Mac OSX

```
git clone https://github.com/ZWang-Lab/LBRAT.git

R CMD INSTALL LBRAT

```
Alternatively, one can install it in R using following code.
### Install LBRAT in R
```
library(devtools)
devtools::install_github(repo = 'ZWang-Lab/LBRAT')

```

## Usage instructions

All functions and examples in the LBRAT are available in the manual (https://github.com/ZWang-Lab/LBRAT/LBRAT-manual.pdf).
