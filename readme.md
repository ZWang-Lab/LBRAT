# LBRAT: R package for retrospective genetic association analysis of longitudinal binary outcomes in GWAS

LBRAT is an R program that performs retrospective association testing for longitudinal binary traits in population samples. 


## Installation Instructions:

### Required software and packages

1. R (http://www.r-project.org/)
  
2. Package nleqslv, lme4, geepack, mvtnorm, data.table, plyr
  
3. PLINK 1.0 or 1.9 (https://www.cog-genomics.org/plink2)

Install the required R package before you install LBRAT package. Install the **LBRAT** using the following steps.


### Install LBRAT on LUNIX or Mac OSX

```
git clone https://github.com/ZWang-Lab/LBRAT.git

R CMD INSTALL LBRAT

```
Alternatively, one can install it in R using the following code.
### Install LBRAT in R
```
library(devtools)
devtools::install_github(repo = 'ZWang-Lab/LBRAT')

```

## Usage instructions:

Details for functions and examples in the LBRAT are available in the manual (https://github.com/ZWang-Lab/LBRAT/blob/master/LBRAT-manual.pdf) and in the vignettes(https://zwang-lab.github.io/LBRAT/articles/L-BRAT.html).

## Citations:

Wu, W., Wang, Z., Xu, K., Zhang, X., Amei, A., Gelernter, J., ... & Wang, Z. (2019). Retrospective Association Analysis of Longitudinal Binary Traits Identifies Important Loci and Pathways in Cocaine Use. Genetics, 213(4), 1225-1236.(https://doi.org/10.1101/628180)
