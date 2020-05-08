####################################################################################
#     Modelling fish habitat selection using functional regression 🐟
#     Created by :    Jeremie Boudreault (Jeremie.Boudreault@ete.inrs.ca)
#     Current file :  Data_salmon_lengths.R
#     Aim :           Convert salmon lengths into number of fish (0+, 1+, 2+)
####################################################################################

#### Part 0 : Libraries, functions and variable definition

# Import library
library('dplyr') # Data.frame manipulation
library('ggplot2') # Data visualisation

# Import raw data of SMR
SMR.raw <- readRDS('data/SMR_2017_field_data_clean.Rda')

###### Part 1 : Lengths of salmon visualisation

# Extract the length from the data.frame
SMR.lengths <- as.vector(na.omit(unlist(SMR.raw[, paste0('Length', 1:11)])))

# Create a data.frame
All.lengths <- data.frame(River = rep('SMR', length(SMR.lengths)), Length = SMR.lengths)

# Producing the graph of the length for the SMR
print(ggplot(data=All.lengths, aes(x=Length)) + 
        geom_histogram(binwidth = 2.5, boundary=0) +
        xlab('Juvenile salmons lengths (mm)') +
        scale_x_continuous(breaks=2:14 * 10) +
        ggtitle('Juvenile salmons lengths distribution'))

# Extract minimum and maximum
All.lengths %>% group_by(River) %>%
                summarize(min = min(Length), max = max(Length))

# Limits of the SMR
Lengths.lim <- data.frame(River = c('SMR'), FryMin = c(30), Parr1Min = c(55), Parr2Min = c(88), Parr2Max = c(125))

# Producing the graph of the lengths + Limits
pdf('out/data visualisation/Figure_2_Salmon_lengths.pdf', width=6, height=4)
print(ggplot(data=All.lengths, aes(x=Length)) +
  geom_histogram(aes(group=River), boundary=0, binwidth=2.5,  fill='#c1c1c1', color='#636363', lwd=0.25) +
  xlab('Length (mm)') +
  ylab('Number of juvenile salmon') +
  scale_x_continuous(breaks=3:14 * 10) +
  scale_y_continuous(limits = c(0, 85)) +
  geom_vline(data = data.frame(x=unlist(Lengths.lim[, 3:5])), aes(xintercept=x), color='red', lwd=1, lty=3) +
  geom_text(data = data.frame(x=c(40, 70, 105), y=c(80, 35, 15), lab=c('Fry (0+)', 'Parr (1+)', 'Parr (2+)')), aes(x=x, y=y, label=lab), size=3)
)
dev.off()

###### Part 2 : Create a new variable to convert lengths to number of salmon

# Merge the lengths limit and the data (SMR)
SMR.clean <- SMR.raw %>%
        mutate(River = 'SMR') %>%
        left_join(Lengths.lim, by='River') 

# Calculate number of Fry, Parr 1+ and Parr 2+ from lengths at SMR
SMR.clean$nFry_M <- sapply(1:nrow(SMR.clean), function(w) sum(SMR.clean[w, paste0('Length', 1:11)] < SMR.clean$Parr1Min[w], na.rm=T))
SMR.clean$nParr1_M <- sapply(1:nrow(SMR.clean), function(w) sum(SMR.clean[w, paste0('Length', 1:11)] >= SMR.clean$Parr1Min[w] & SMR.clean[w, paste0('Length', 1:11)] < SMR.clean$Parr2Min[w], na.rm=T))
SMR.clean$nParr2_M <- sapply(1:nrow(SMR.clean), function(w) sum(SMR.clean[w, paste0('Length', 1:11)] >= SMR.clean$Parr2Min[w] & SMR.clean[w, paste0('Length', 1:11)] < SMR.clean$Parr2Max[w], na.rm=T))
SMR.clean$nParr3_M <- sapply(1:nrow(SMR.clean), function(w) sum(SMR.clean[w, paste0('Length', 1:11)] >= SMR.clean$Parr2Max[w], na.rm=T))

# Validation
SMR.clean %>% group_by(River) %>% summarize(nFry = sum(nFry_M), 
                                          nParr1 = sum(nParr1_M), 
                                          nParr2 = sum(nParr2_M),
                                          nParr3 = sum(nParr3_M))
  
# Calculate the total number of fish
SMR.clean <- SMR.clean %>% 
  mutate(nFry_tot = as.integer(nFry + nFry_M),
         nParr_M = as.integer(nParr1_M + nParr2_M + nParr3_M),
         nParr1_tot = as.integer(nParr1_M),
         nParr2_tot = as.integer(nParr2_M),
         nParr_tot = as.integer(nParr + nParr_M),
         nSalmon_M = as.integer(nFry_M + nParr_M), 
         nSalmon_tot = as.integer(nFry_tot + nParr_tot + nSalmon))

# Spot check on the summry
SMR.clean %>% group_by(River) %>% summarize(nFry = sum(nFry), 
                                          nFry_M = sum(nFry_M),
                                          nFry_tot = sum(nFry_tot),
                                          nParr1_M = sum(nParr1_M), 
                                          nParr1_tot = sum(nParr1_tot),
                                          nParr2_M = sum(nParr2_M),
                                          nParr2_tot = sum(nParr2_tot),
                                          nParr3_M = sum(nParr3_M),
                                          nParr_M = sum(nParr_M),
                                          nParr = sum(nParr),
                                          nParr_tot = sum(nParr_tot),
                                          nSalmon_M = sum(nSalmon_M),
                                          nSalmon = sum(nSalmon),
                                          nSalmon_tot = sum(nSalmon_tot))

# Here, we will only keep the variable of interest and reorganise the data.frame
names(SMR.clean)

# Cleaning of the dataset
SMR.clean <- SMR.clean %>% 
             dplyr::select(Site, Parcelle, nFry_M, nFry_tot, nParr1_M, nParr1_tot, nParr2_M, nParr2_tot, nParr_M, nParr_tot, nSalmon_M, nSalmon_tot, Velocity, Depth, Temp, D50, D84)

# Verify is there is any NAs
which(is.na(SMR.clean), arr.ind=T)

# Verify the class of the variables
sapply(SMR.clean, class)

# Save the dataset
saveRDS(SMR.clean, 'data/SMR_2017_field_data_clean_nsalmon.Rda')
