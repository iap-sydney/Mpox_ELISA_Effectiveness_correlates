# Estimating and predicting effectiveness of vaccinia vaccination against Mpox: A meta-analysis and modelling study

This repository contains code to recreate all analysis presented in the corresponding article (hyperlink to be added)

## Cloning repository and setting up RStudio

This project can be cloned from github into RStudio (or other IDE) using instructions [here](https://happygitwithr.com/existing-github-first.html).

Once cloned I suggest changing the project settings in RStudio under `Tools` \> `Project Settings` :

-   Restore .RData workspace at startup: `No`

-   Save workspace to .RData on exit: `No`

This ensures you're starting from a clean session whenever you open the project in RStudio. You can set similar settings in `Global options` which will affect all projects.

## Running the project

The project is to run from start to finish by executing the `driver.R` script. This script should always run cleanly (without error) on the `main` branch.

`driver.R` calls other scripts including `setup.R` and files in the `processing/` and `analysis/` folders.

All data used during the project is contained in the `data/` folder

## Packages required 
-   rstan
-   ggpubr
-   tidyverse
-   loo
-   latex2exp
-   nleqslv
