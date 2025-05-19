# Refining Relational Event Models: Bayesian Penalization and Variable Selection in REMs

## Introduction

This repository contains the data and R scripts to reproduce the findings for the master's thesis *Refining Relational Event Models: Bayesian Penalization and Variable Selection in REMs*, which investigates the use of Bayesian regularization techniques for variable selection and predictions in Relational Event Models (REMs).

REMs are powerful tools for modeling dynamic interaction data over time, commonly used in the analysis of relational event history (REH) data. While their flexibility is an advantage, REMs face challenges in variable selection due to the abundance of potentially collinear endogenous and exogenous predictors. This complexity risks overfitting and reduces interpretability, especially in high-dimensional settings.

To address these issues, this study systematically compares Exact Bayesian Regularization (EBR) through a reparameterization of REMs to Poisson regressions and Approximate Bayesian Regularization (ABR) normally approximating the Likelihood function using Ridge and Horseshoe priors against standard Maximum Likelihood Estimation (MLE). Through simulation studies and an empirical application to Spotify collaboration data, we evaluate these methods in terms of variable selection, bias, variance, and predictive performance. The results show that ABR with suitable priors provides a robust and computationally efficient alternative to EBR and outperforms MLE, particularly for small sample sizes.

## Data

This repository contains data in the following structure:

``` plaintext
├── BayesianRegularization-REMs.Rproj
├── Data
│   ├── 01_simulation
│   └── 02_application
├── Output
│   ├── 01_result_files
│   ├── 02_plots
│   └── 03_tables
├── README.md
├── Script
│   ├── 01_functions
│   └── 02_analyses
├── renv
│   ├── activate.R
│   └── library
├── renv.lock
└── Requirements.md
```

### `BayesianRegularization-REMs.Rproj`

This is the R project file for the repository. It contains the project settings and configuration. When running the project, it is useful to open this file in RStudio.

### `Data`

This folder contains the data used in the analyses. It is divided into two subfolders.

#### `01_simulation`

This folder contains the simulated data used in the analyses.

We generated **100** directed REH datasets, each containing **7,400 events** among **50 actors**. For each dataset, we created subsets by truncating to the first **100**, **200**, **400**, **800**, **1,600**, **3,200**, and **6,400** events, leading to **700** datasets in total for the main simulation. These datasets can be found in `Data/01_simulation/01_edgelists_dependent`.

To assess robustness, we additionally generated independent datasets of each size **M + 1000**. These can be found in `Data/01_simulation/01_edgelists_independent`.

The datasets were generated using the following effects:

-   **Endogenous effects**: 3 non-zero effects (*reciprocity*, *indegreeSender*, and *outdegreeReceiver*) from the `remstats` package.
-   **Exogenous effects**: 12 non-zero *Min* and *Max* effects from 6 exogenous covariates (3 continuous and 3 binary) with effects of varying strengths

Using the `remulate` package, dyads were sampled iteratively from a fixed risk set of 2,450 possible possible sender-receiver pairs.

The directory additionally includes the following files:

-   `covar.RData`: This file contains a data frame with the generated exogenous covariates of the 50 actors.
-   `parameters.RData`: This file contains the list of data-generating parameters.

#### `02_application`

The folder `02_application` contains the data used in the empirical application. The data is stored in `.RData` format and includes:

-   `Data_Popular_100.RData`: This file contains dates and actors of collaborations between the most popular 100 Spotify artists that collaborated with other artists at least twice from 2010 to 2023. There are a total of **244** collaboration events.

-   `covar_spotify.RData`: This file contains a data frame with the exogenous covariates of the 100 artists. These include gender, country of origin, popularity, and age.

### `Output`

The `Output` folder contains the results of the analyses. It is divided into three subfolders.

#### `01_result_files`

This folder contains the result files of the analyses. The results are stored in `.RData` and include the estimates of the models, the selected variables following the explored selection criteria, the discovery rates, the distance metrics, Matthews' correlation coefficient (MCC), and the predictive performance metrics. Since the results are not nicely formatted, they are summarized into plots and tables in the next two folders.

#### `02_plots`

This folder contains the plots of the results. The plots are stored in `.png` format and stored in the following subfolders:

-   `01_coefficient_plots`: This folder contains plots comparing the extent of bias and variance of the estimates from ABR and MLE from both the dependent and independent data. In the plots, the estimates are shown as black dots, and the true data-generating parameters are shown as red dots. The 95% credible intervals (CIs) for the estimates are displayed as ribbons.

-   `02_variable_selection`: This folder contains plots comparing the selection performance of the explored selection criteria. Plots include the True Discovery Rate (TDR), False Discovery Rate (FDR), the distance metric, and Matthews' correlation coefficient (MCC) for the dependent and independent data.

-   `03_predictive_performance`: This folder contains plots comparing the predictive performance of full models (including MLE, ABR and EBR) and sparse models after applying the explored selection criteria. Plots are created for both in-sample (everything starting with `pp_is`) and out-of-sample (everything starting with `pp_oos`) predictive performance metrics.

#### `03_tables`

This folder contains the tables of the results. The tables are stored in `.csv` format and include the bias and variance of the estimates from the models, the discovery rates, the distance metrics, Matthews' correlation coefficient (MCC), and the predictive performance metrics. The tables are in the `.tex` format.

### `Script`

The `Script` folder contains the R scripts used in the analyses. It is divided into two subfolders: `01_functions`, where the functions are stored, and `02_analyses`, which contains the scripts for the analyses sourcing the functions from `01_functions`. Further information about the Script can be found in *Reproducing Results* below.

### `renv` and `renv.lock`

The `renv` folder and the `renv.lock` file contain information about the R environment used in the analyses. Through the renv package, the R environment can be restored to the same state as when the analyses were run.

### `Requirements.md`

The `Requirements.md` file contains the software requirements and refers to the `renv` file that contains the package versions used in the analyses.

## Reproducing Results

### Prerequisites

To reproduce the analyses, the following are needed:

-   **R** (≥ 4.4.2)
-   The **`renv`** package for managing the R environment
-   **CmdStanR** backend for Bayesian model fitting
-   Other software and dependencies listed in `Requirements.md`

> ⚠️ **Hardware Requirements:**\
> Analyses were run on a high-performance machine with the following specifications:
>
> -   **CPU**: 112 × Intel(R) Xeon(R) Platinum 8580
> -   **Threads**: 224
> -   **Memory**: 851 GB RAM
>
> Due to processing large arrays containing endogenous statistics for all potential dyads in the risk set for all time points, particularly this memory is needed in order to allow parallel processing.

Specific package versions can be found in `Requirements.md`.

### Running the Script

The results can be reproduced by running the script as outlined in the steps below:

1.  Open the R project file `BayesianRegularization-REMs.Rproj` in RStudio.

2.  Install the required packages by running the following command in the R console:

    ``` r
    renv::restore()
    ```

    Additionally, install the `cmdstanR` backend by running the following command in the R console:

    ``` r
    cmdstanr::install_cmdstan()
    ```

3.  Open and run the R script `Script/02_analyses/01_generating_data.qmd`. This script generates the data frame with covariates `covar.RData` and the list of data-generating parameters `parameters.RData` for the simulation study. Following this, using the `remulate` and `future` packages, the script generates the datasets for the simulation study. The datasets are saved in `Data/01_simulation/01_edgelists_dependent` and `Data/01_simulation/01_edgelists_independent`. Lastly bias and variance of the estimates from the MLE models are evaluated and saved in tables in `Output/03_tables` and plots in `Output/02_plots/01_coefficient_plots`.\
    **Note:** Depending on the number of available cores, the script will take considerable time to run. With the specifications mentioned above, it took approximately **5 hours** to generate the data.

4.  Run the script in `Script/02_analyses/02a_estimating_models_dependent.qmd` and `Script/02_analyses/02b_estimating_models_independent.qmd` to estimate the MLE, ABR and EBR models on the dependent and independent datasets, respectively. The results are saved in `Output/01_result_files/01a_estimates_dependent`, `Output/01b_result_files/01b_estimates_ebr` and `Output/01_result_files/01c_estimates_independent`.\
    **Note:** The estimation of the models, particularly that of the EBR models using `brms`, will take considerable time to run. With the specifications mentioned above, it took multiple days to estimate all models.

5.  Run the script in `Script/02_analyses/03a_selecting_variables_dependent.qmd` and `Script/02_analyses/03b_selecting_variables_independent.qmd` to select the variables using the explored selection criteria and explore that selection through the computation of the discovery rates, distance metrics and Matthews' correlation coefficient (MCC). To that end, the script also generates plots and tables for the results. The results are saved in `Output/01_result_files/02a_selection_dependent`, `Output/01_result_files/02b_selection_independent`, `Output/01_result_files/03a_discovery_rates_dependent`, `Output/01_result_files/03b_discovery_rates_independent`, and `Output/02_plots/02_variable_selection` and `Output/03_tables`.\
    **Note:** The selection of variables, i.e., the first code chunk takes moderate time to run (\~20-30 minutes). After the selection, the script will be fast to run.

6.  Run the script in `Script/02_analyses/04a_predicting_performance_dependent.qmd` and `Script/02_analyses/04b_predicting_performance_independent.qmd` to evaluate the predictive performance of the models. The results are saved in `Output/01_result_files/04a_predictive_performance_is_dependent`, `Output/01_result_files/04b_predictive_performance_oos_dependent`, `Output/01_result_files/04c_predictive_performance_is_independent`, `Output/01_result_files/04d_predictive_performance_oos_independent`, `Output/02_plots/03_predictive_performance` and `Output/03_tables`.\
    **Note:** The computation of the predictive performance metrics takes multiple hours to run. After the computation, the plots and tables are generated quickly.

7.  Run the script in `Script/02_analyses/05_application.qmd` to estimate the models on the Spotify collaboration data, illustrate variable selection and evaluate the predictive performance. The resulting estimates are saved in `Output/01_result_files/05_application` and the plot in `Output/02_plots/04_application`. 
  **Note:** The estimation of the models using MLE and ABR is fast, however, the estimation of the EBR models using `brms` will take considerable time to run (approximately 5 hours).
  
## Ethics and Privacy

Ethics approval was granted by [Ethics Review Board of the Faculty of Social & Behavioural Sciences at Utrecht University](https://ferb.sites.uu.nl). The ethical approval case numbers are 24-2053, 24-2054, and 24-2057.

Due to the synthetic nature of the data used in the simulation study, no privacy concerns are relevant. The simulated data is generated from a known model and does not contain any real-world personal data.

By using the Spotify collaboration data, we acknowledge that the data is publicly available and does not contain any personally identifiable information. The data was collected from the Spotify API and is used in accordance with Spotify's terms of service.

## License

This project is licensed under the GNU General Public License v3.0. See the [LICENSE](LICENSE.md) file for details.

## Permissions and Access

This archive will indefinitely be publicly available on [GitHub](https://github.com/jonathankoop/BayesianRegularization-REMs). Full responsibility for the content of this archive lies with [Jonathan Koop](https://jonathankoop.eu/). In the case of questions, do not hesitate to contact me by emailing [j.koop@uu.nl](mailto:j.koop@uu.nl) or [jonathankoop@proton.me](mailto:jonathankoop@proton.me).
