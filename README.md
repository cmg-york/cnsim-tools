![CNSim-Bitcoin Logo](src/site/resources/images/logo.png)

# Data Analysis Tools for CNSim (cnsim-tools).

This repository contains tools for analyzing [CNSim](https://github.com/cmg-york/cnsim-engine/) logs. Python-based tools are under `Py-Tools` and R-based tools are under `R-Tools`.

Individual tools are under construction or deprecation, except for those documented below.

## R-Tools
(not to be confused with [RTools](https://cran.r-project.org/bin/windows/Rtools/) for building R packages)

### logAnalysis

A set of tools for belief and pace analysis of CNSim logs.
- `library.R`: contains all functions for performing belief and pace analysis of bitcoin data obtained via e.g. [cnsim-bitcoin](https://github.com/cmg-york/cnsim-engine/). Include in your R script via `source('library.R')`.
- Refer to the comments on top of each function in `library.R` for documentation.