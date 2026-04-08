# The combined effects of ultraviolet-B radiation and temperature on the survival of embryos from three tropical anuran species

This repository contains data and R scripts used to analyze the combined effects of ultraviolet-B radiation (UVBR) and temperature on embryo survival and malformations in three tropical anuran species.

## Author

Katalina Gutiérrez Hernández  
Universidad del Tolima

## Repository contents

Scripts/  
&nbsp;&nbsp;stats/  
&nbsp;&nbsp;&nbsp;&nbsp;mortality_analysis.R  
&nbsp;&nbsp;&nbsp;&nbsp;malformations_analysis.R  
&nbsp;&nbsp;figures/  
&nbsp;&nbsp;&nbsp;&nbsp;mortality_figures.R  
&nbsp;&nbsp;&nbsp;&nbsp;malformations_figures.R  
Data/  
&nbsp;&nbsp;mortality_counts.csv  
&nbsp;&nbsp;malformations_counts.csv  

## Study overview

The analyses in this repository evaluate the effects of UVBR and temperature treatments on embryos of three tropical anuran species:

- *Boana platanera*
- *Engystomops pustulosus*
- *Rhinella horribilis*

Two response variables were analyzed:

- embryo survival
- embryo malformations

## Scripts

### Statistical analyses

`Scripts/stats/mortality_analysis.R`  
Statistical analyses for embryo survival.

`Scripts/stats/malformations_analysis.R`  
Statistical analyses for embryo malformations.

### Figures

`Scripts/figures/mortality_figures.R`  
Scripts for observed and predicted survival figures.

`Scripts/figures/malformations_figures.R`  
Scripts for observed and predicted malformation figures.

## Data

`Data/mortality_counts.csv`  
Data used for survival analyses.

`Data/malformations_counts.csv`  
Data used for malformation analyses.

## Requirements

The scripts were written in R and use the following packages:

- `dplyr`
- `ggplot2`
- `lme4`
- `emmeans`

## Notes

- Survival and malformation analyses were performed separately.
- Inferential models excluded the UVBR = 0 treatment when required to avoid complete separation.
- Figure scripts include options for color and black-and-white versions.
