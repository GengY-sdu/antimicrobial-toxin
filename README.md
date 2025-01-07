# Global biogeography and projection of antimicrobial toxin genes
## Introduction

This project includes all R codes for machine learning analysis of antimicrobial toxin genes (ATGs).

Please place these files in a folder named `antimicrobial_toxin` and set default working directory of R to the previous level of the `antimicrobial_toxin`. For example, if the path of `antimicrobial_toxin` is `a/b/c/antimicrobial_toxin`, execute the command `setwd(dir = "a/b/c/")` in R.

## Detailed description

R codes involve four processes of machine learning: feature selection, hyperparameter tuning, test model, and modelling.

-   Folder `rf_cluster`: R codes used for predicting the global distribution of `ATG cluster diversity`.

-   Folder `rf_family`: R codes used for predicting the global distribution of `ATG family diversity`.

-   Folder `rf_abundance`: R codes used for predicting the global distribution of `ATG abundance`.

-   `train_new.csv`: The dataset used for machine learning.

## Please note

Please download the map data for analysis from **http://bioinfo.qd.sdu.edu.cn/kegg/map_data/map_toxin.zip** using a browser.
