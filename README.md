## Modelling fish physico-thermal habitat selection using functional regression 🐟

* A scientific research by **Jeremie Boudreault**, Fateh Chebana, André St-Hilaire and Normand Bergeron
* This project is part of my **master degree** in water sciences at [**Institut National de la Recherche Scientifique**](http://inrs.ca) (INRS)
* All **codes** and **data** are made available here under the [**Creative Common License** ![](https://i.creativecommons.org/l/by-nc-nd/4.0/80x15.png)](http://creativecommons.org/licenses/by-nc-nd/4.0/) 
* Questions regarding the code or the data should be sent to **Jeremie.Boudreault@ete.inrs.ca**

## Data

Data are from **field survey** that have been conducted during summer 2017 on the  **Sainte-Marguerite river** (SMR) :

* `data/field/*` : contains the raw .xlsx file filled after each day of field work 
* `data/*` : contains the cleaned and transformed datasets 

## R codes

Codes are all from Jeremie Boudreault. They use the R package `mgcv` to fit **generalized additive models** (GAM) and of `FDboost` to fit **functional regression models** (FRM) :

* `R/Data_initial_cleaning.R` : code to clean the field data spreadsheets and produce more adapted datasets
* `R/Data_salmons_lengths.R`: code to convert the salmon lengths to number of fry and parr
* `R/Data_per_site.R` : code to produce the observations at each site (mean value or functional observations)
* `R/GAMs_all.R` : code to fit several types of GLM/GAM on the data using the `mgcv` package
* `R/GAMs_best.R` : among all models, do variable selection to find the best GAM models and save them to `out/models`
* `R/GAMs_predictions` : code to calculate the leave-one-out predictions for the GAMs
* `R/FRMs_all.R` : code to fit several types of FRM on the data using the `FDboost` package
* `R/FRMs_best.R` : among all models, select the best FRM models and save them to `out/models` 
* `R/FRMs_predictions` : code to calculate the leave-one-out predictions for the FRMs
* `R/Models_coefficients.R` : code to extract the coefficients of the best models
* `R/Models_performance.R` : compare the results between GAMs and FRMs

## Out

A folder for the **results** at each part of the coding process :

* `out/data visualisation/*` : raw data visualisation and tables
* `out/models/*` : fitted final models
* `out/coefficients/*` : coefficients of the models
* `out/predictions/*` : leave-one-out predictions and goodness-of-fit
