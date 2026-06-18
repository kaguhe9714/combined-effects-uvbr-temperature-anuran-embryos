# The combined effects of ultraviolet-B radiation and temperature on the survival of embryos from three tropical anuran species

This repository contains data and R scripts used to analyze the combined effects of ultraviolet-B radiation (UVBR) and temperature on embryo survival and to summarize embryo malformations in three tropical anuran species.

## Author

Katalina Gutiérrez Hernández  
Universidad del Tolima

## Repository contents

```text
Scripts/
  stats/
    mortality_analysis.R
  exploratory/
    malformations_exploratory_not_reported.R
  figures/
    mortality_figures.R
    malformations_figures.R

Data/
  mortality_counts.csv
  malformations_counts.csv

Note: The malformations script is provided as an exploratory analysis only. This model was not used for manuscript inference because the number of surviving embryos available for malformation assessment was highly unbalanced among treatments/species. Malformations are reported descriptively in the manuscript.

Study overview

The analyses in this repository evaluate the effects of UVBR and temperature treatments on embryos of three tropical anuran species:

Boana platanera
Engystomops pustulosus
Rhinella horribilis

Two response variables were evaluated:

embryo survival, analyzed inferentially
embryo malformations, summarized descriptively
Scripts
Statistical analyses

Scripts/stats/mortality_analysis.R
Statistical analyses for embryo survival.

Exploratory analyses

Scripts/exploratory/malformations_exploratory_not_reported.R
Exploratory analysis for embryo malformations. This model was not used for manuscript inference because the number of surviving embryos available for malformation assessment was highly unbalanced among treatments/species. Malformations are reported descriptively in the manuscript.

Figures

Scripts/figures/mortality_figures.R
Scripts for observed and predicted survival figures.

Scripts/figures/malformations_figures.R
Scripts for descriptive malformation figures.

Data

Data/mortality_counts.csv
Data used for survival analyses.

Data/malformations_counts.csv
Data used to summarize malformations descriptively and to generate the malformation figure.

Requirements

The scripts were written in R and use the following packages:

dplyr
ggplot2
lme4
emmeans
Notes
Survival analyses and descriptive malformation summaries were performed separately.
The inferential survival model excluded the UVBR = 0 treatment to avoid complete separation.
The malformation model is provided only as an exploratory analysis and was not used for manuscript inference.
Figure scripts include options for color and black-and-white versions.
