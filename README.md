# COVID19-HHH4-Africa

This is the R code associated with the paper "Tracking and predicting the African COVID-19 pandemic" ([medrxiv link](https://www.medrxiv.org/content/10.1101/2020.11.13.20231241v1))

This code was tested with R version 3.6.3 on Ubuntu 18.04 LTS and Windows and version 4.0.2 on Mac OS

Set the working directory to the one with these source files, and run these files:

- `01_data-processing-model.R`: cleans and processes the input data in `data/original` and gets it ready for modelling. The output of this script goes into the `data/processed` folder.

- `01_data-processing-plot.R`: cleans and processes the input data in `data/original` and gets it ready for plotting. The output of this script goes into the `data/processed` folder.

- `02_modelfitting.R`: uses the processed data to fit the models and saves them in `output/models/`. Summary tables are also produced and saved in `output/tables`.

- `03_figures-model.R`: uses fitted model to create figures in the `figs` folder.

- `03_figures_aux.R`: use processed data and fitted model to create other figures in the `figs` folder.


