# Viral disease outcomes are indistinguishable between experimentally infected bats and rodents

Data and code to reproduce the analyses in ["Viral disease outcomes are indistinguishable between
experimentally infected bats and rodents"](https://doi.org/10.1101/2025.08.28.672855)

## Package Management

Run `renv:restore()` (see `packages.R`) in the root directory to load the correct package versions.

Then set working directory to `scripts` and run scripts in the following order:


## Data Preprocessing

1. `virus_taxize.R` - grabs virus higher level taxonomy from NCBI 
2. `virus_traits_merge.R` - merges reservoir and fatality data from Mollentze and Guth papers
3. `host_specificity.R` - gets host range data from VIRION and calculates host specificity metrics and host evolutionary isolation
4. `symptoms_harmonize.R` - cleans raw data on organs damaged, behavioural changes, and symptoms. Creates severity scores. Note that species only and individual level data are merged in the output (Individual level data have values in the IndividualID column; species only data have NA). Makes summary plots for individuals.
5. `dose_harmonize.R` - cleans dose volume, dose unit, and inocculation route for individual data. Marges with body size data from COMBINE to calculate dose/mass. Makes summary plots.  


## Data Summaries 

6. `data_summaries.R` - calculates simple quantiateive measures for reporting in the paper. Creates heatmap of host x virus severity.
7. `phylo_plot_sampling.R` - generates descriptive figures for phylogenetic sampling, host traits, and numbers of described and experimental hosts per virus. 


## Analyses

8. `brms_models.rmd` - conducts Bayesian hierarchical models via brms + Stan and outputs `brms_models.html`

9. `model_summaries.rmd` - generates combined forest plots for brms models. 

10. `scaling_plotting_functions.R` - common functions for formatting data and outputs (plots and tables).

## Data

`raw_data` contains untransformed original data from this manuscript (`article_metadata.csv`, `primary_screen.csv`, `species_only_data.csv`, `individual_data.csv`), taxonomic name translation tables, and published data (e.g. [VIRION](https://journals.asm.org/doi/10.1128/mbio.02985-21), [Guth et al. (2022)](https://www.pnas.org/doi/10.1073/pnas.2113628119), [Mollentze & Streicker (2020)](https://www.pnas.org/doi/full/10.1073/pnas.1919176117), [Upham et al. 2019](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000494), [COMBINE](https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecy.3344)), and Web of Science search results and relevant [CLOVER](https://academic.oup.com/bioscience/article/71/11/1148/6353869) articles included in the primaty literature search. 


`clean_data` contains cleaned and summarized versions of the `raw_data` files, including harmonized dose and severity measures, host range data, viral traits, and virus taxonomy.

*Note: output .rds files in `fit_models` are not included to reduce repository size, but will be re-fit and saved when `brms_models.rmd` is compiled.* 