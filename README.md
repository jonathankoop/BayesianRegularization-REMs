# Refining Relational Event Models: Bayesian Penalization and Variable Selection in REMs

## Introduction

This repository contains the data and R scripts to reproduce the findings for the master's thesis *Refining Relational Event Models: Bayesian Penalization and Variable Selection in REMs*, which investigates the use of Bayesian regularization techniques for variable selection and predictions in Relational Event Models (REMs).

REMs are powerful tools for modeling dynamic interaction data over time, commonly used in the analysis of relational event history (REH) data. While their flexibility is an advantage, REMs face challenges in variable selection due to the abundance of potentially collinear endogenous and exogenous predictors. This complexity risks overfitting and reduces interpretability, especially in high-dimensional settings.

To address these issues, this study systematically compares Exact Bayesian Regularization (EBR) through a reparameterization of REMs to Poisson regressions and Approximate Bayesian Regularization (ABR) normally approximating the Likelihood function using Ridge and Horseshoe priors against standard Maximum Likelihood Estimation (MLE). Through simulation studies and an empirical application to Spotify collaboration data, we evaluate these methods in terms of variable selection, bias, variance, and predictive performance. The results show that ABR with suitable priors provides a robust and computationally efficient alternative to EBR and outperforms MLE, particularly for small sample sizes.

## Data

This repository provides both simulated and empirical data:

### Simulated Data

We generated **100** directed REH datasets, each containing **7,400 events** among **50 actors**. For each dataset, we created subsets by truncating to the first **100**, **200**, **400**, **800**, **1,600**, **3,200**, and **6,400** events, leading to **700** datasets in total for the main simulation. These datasets can be found in `Data/01_simulation/01_edgelists_dependent`.

To assess robustness, we additionally generated independent datasets of each size **M + 1000**. These can be found in `Data/01_simulation/01_edgelists_independent`.

For out-of-sample predictive evaluation, the following **1,000 events** of each dataset were used.

- **Endogenous effects**: 3 non-zero effects (*reciprocity*, *indegreeSender*, and *outdegreeReceiver*) from the `remstats` package.
- **Exogenous effects**: 12 non-zero *Min* and *Max* effects from 6 exogenous covariates (3 continuous and 3 binary) with effects of varying strengths

Using the `remulate` package, dyads were sampled iteratively from a fixed risk set of 2,450 possible possible sender-receiver pairs.

### Empirical Data

The empirical example uses collaborations between popular Spotify artists from 2010 to 2023. This dataset includes:

- **244 directed collaboration events**
- **62 unique artists**
- **Exogenous attributes**: Gender, Age, Country, and Popularity (from Spotify and Last.fm)
- **Endogenous statistics**: All directed REM statistics provided by the `remstats` package

## Reproducing Results

### Prerequisites

To reproduce the analyses, the following are needed:

- **R** (≥ 4.4.2)
- The **`renv`** package for managing the R environment
- **CmdStanR** backend for Bayesian model fitting
- Other software and dependencies listed in `Requirements.md`

> ⚠️ **Hardware Requirements:**  
> Analyses were run on a high-performance machine with the following specifications:
>
> - **CPU**: 112 × Intel(R) Xeon(R) Platinum 8580  
> - **Threads**: 224  
> - **Memory**: 851 GB RAM  
>
> Due to processing large arrays containing endogenous statistics for all potential dyads in the risk set for all time points, particularly this memory is needed in order to allow parallel processing.

### Running the Script

The results can be reproduced by running the script as outlined in the steps below:

1. 
