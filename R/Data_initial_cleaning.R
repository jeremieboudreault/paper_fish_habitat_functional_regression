####################################################################################
#     Modelling fish habitat selection using functional regression 🐟
#     Created by :     Jeremie Boudreault (Jeremie.Boudreault@ete.inrs.ca)
#     Current file :   Data_initial_cleaning.R
#     Aim :            Clean the field data and produce nice files to use
####################################################################################

##### Part 0 : Libraries, functions and variable definition

# Import library
library('xlsx') # Importing Excel files (Field data)
library('dplyr') # Working with data.table easily

##### Part 1 : Data importation

# Importation of SMR data
SMR.raw <- read.xlsx('data/Field/Field_data_summer_2017_SMR.xlsx', 1, startRow = 3, header = T)

# Check the class of the variables
sapply(SMR.raw, class)

# Some transformation of the class
SMR.raw <- SMR.raw %>%
  filter(SiteNew != 999) %>%       # Remove a site where field work was imprecise
  mutate(Site = as.integer(SiteNew),
         Transect = as.integer(Transect),
         Parcelle = as.integer(Parcelle),
         GPS = as.integer(GPS),
         areTD = as.integer(areTD),
         areMB = as.integer(areMB),
         nFry = as.integer(nFry),
         nParr = as.integer(nParr),
         nSalmon = as.integer(nSalmon),
         nRHCA = as.integer(nRHCA),
         nSAFO = as.integer(nSAFO),
         nUnknown = as.integer(nUnknown)) %>%
  dplyr::select(-SiteNew) %>%
  arrange(Site, Parcelle)

# Replaces NAs by 0 in number of fishes
SMR.raw$nFry[is.na(SMR.raw$nFry)] <- 0
SMR.raw$nParr[is.na(SMR.raw$nParr)] <- 0
SMR.raw$nSalmon[is.na(SMR.raw$nSalmon)] <- 0

# Save the file in the data folder
saveRDS(SMR.raw, 'data/SMR_2017_field_data_clean.Rda')
