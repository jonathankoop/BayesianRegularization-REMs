# Software Requirements

## R Version

- R version 4.4.2 (2024-10-31) - "Pile of Leaves"

## R Packages

A full list of all R packages used in the project can be found in the [`renv.lock`](renv.lock) file. The following packages are required to be installed to run the code in this repository:

- `renv` version 1.1.4

To install the required packages, install the `renv` version above and run the following command in R:

```r
renv::restore()
```

This then loads all packages listed in the [`renv.lock`](renv.lock) file into the R environment.