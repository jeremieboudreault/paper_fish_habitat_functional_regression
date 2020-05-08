####################################################################################
#     Modelling fish habitat selection using functional regression 🐟
#     Created by :    Jeremie Boudreault (Jeremie.Boudreault@ete.inrs.ca)
#     Current file :  Data_per_site.R
#     Aim :           Combine data to have functional observations and mean values per site
####################################################################################

#### Part 0 : Libraries, functions and variable definition

# Import library
library('tidyr') # To pivot data.frame
library('dplyr') # Data.frame manipulation
library('ggplot2') # Data visualisation
library('ggpubr') # To have multiple plot on the same graph
library('reshape') # For data.frame melting
library('fitdistrplus') # For fitting distributions to data

# Import clean data of SMR
SMR.clean <- readRDS('data/SMR_2017_field_data_clean_nsalmon.Rda')

# Variable - nKnots :: the number of knots for the density estimate
nKnots <- 2^9

# Variable - Vars :: the variable we are dealing with
Vars <- c('Velocity', 'Depth', 'D50', 'D84', 'Temp')

# Function - findRange :: to find the range of a vector c(xmin, xmax)
findRange <- function(x, liminf="yes", factor=1.1) {
  if (liminf == "yes")
    liminf = min(x, na.rm=T)/factor
  return(c(liminf, max(x, na.rm=TRUE) * factor))
}

#### Part 1 : Generate dataset per site

# Convert site to integer
SMR.clean$Site <- as.integer(SMR.clean$Site)

# Generate the data for the SMR
SMR.per.site <- SMR.clean %>%
                group_by(Site) %>%
                summarize(nFry_M = sum(nFry_M),
                        nFry_tot = sum(nFry_tot),
                        nParr1_M = sum(nParr1_M),
                        nParr1_tot = sum(nParr1_tot),
                        nParr2_M = sum(nParr2_M),
                        nParr2_tot = sum(nParr2_tot),
                        nParr_M = sum(nParr_M),
                        nParr_tot = sum(nParr_tot),
                        nSalmon_M = sum(nSalmon_M),
                        nSalmon_tot = sum(nSalmon_tot),
                        Depth = mean(Depth),
                        Velocity = mean(Velocity),
                        D50 = mean(D50),
                        D84 = mean(D84),
                        Temp = mean(Temp))

# Save the data (SMR)
saveRDS(SMR.per.site, 'data/SMR_2017_clean_per_site.Rda')

#### Part 2 : Generate a general table of information about the site

# Start by melting the data.frame of the complete observation
SMR.clean.var.melt <- melt(SMR.clean, id.vars = c('Site', 'Parcelle'), measure.vars = c('Depth', 'Velocity', 'D50', 'Temp'))

# Melting the data.frame of the per site observation
SMR.per.site.melt <- melt(as.data.frame(SMR.per.site), id.vars = c('Site'), measure.vars = c('Depth', 'Velocity', 'D50', 'Temp'))

# We create a summary of the physical variable
SMR.var.summary <- SMR.clean.var.melt %>% 
  group_by(Site, variable) %>%
  # summarize(Mean = paste0(round(mean(value), 1), ' (±', round(sd(value), 1), ')')) %>%
  summarize(Mean = paste0(round(mean(value), 4))) %>%
  pivot_wider(id_cols = Site, names_from=variable, values_from = Mean) 

# We create a grand summary across site
SMR.per.site.summary <- SMR.per.site.melt %>%
  group_by(variable) %>%
  summarize(Mean = paste0(round(mean(value), 3), ' (±', round(sd(value), 3), ')')) %>% # With standard deviation
  pivot_wider(names_from = variable, values_from=Mean) %>%
  mutate(Site = 0, 
         nFry_M = paste0(round(mean(SMR.per.site$nFry_M), 1), ' (±', round(sd(SMR.per.site$nFry_M), 1), ')'), 
         nParr1_M = paste0(round(mean(SMR.per.site$nParr1_M), 1), ' (±', round(sd(SMR.per.site$nParr1_M), 1), ')'),
         nParr2_M = paste0(round(mean(SMR.per.site$nParr2_M), 1), ' (±', round(sd(SMR.per.site$nParr2_M), 1), ')')) 

# We merge the SMR.var.summary with the number of fish per site
SMR.summary <- SMR.var.summary %>%
  left_join(SMR.per.site[, c('Site', 'nFry_M', 'nParr1_M', 'nParr2_M')], by='Site') %>%
  mutate(nFry_M = as.character(nFry_M), nParr1_M = as.character(nParr1_M), nParr2_M = as.character(nParr2_M)) %>%
  bind_rows(SMR.per.site.summary)   %>%          
  dplyr::rename(Fry = nFry_M, Parr1 = nParr1_M, Parr2 = nParr2_M)

# Save to .csv
write.csv2(SMR.summary, 'out/data visualisation/Table_1_Summary_across_site.csv')

#### Part 3 : Detemine to which family each response variable corresponds

# We will focus on the Fry and Parr1 only
Y.all <- c('nFry_M', 'nParr1_M')

# Create two data.frame to store the densities and the AIC
densities <- data.frame(Y='NA', x=NA, y.norm=NA, y.pois=NA, y.nbinom=NA)
AIC <- data.frame(Y=NA, AIC.norm=NA, AIC.pois=NA, AIC.nbinom=NA)

# Start loop on all Y value
for (Y in Y.all) {
  
  # Fit normal, poison and negative binomial to distributions (first Y)
  fit.norm <- fitdist(unlist(SMR.per.site[, Y]), distr=c('norm'))
  fit.pois <- fitdist(unlist(SMR.per.site[, Y]), distr=c('pois'))
  fit.nbinom <- fitdist(unlist(SMR.per.site[, Y]), distr=c('nbinom'))
  
  # Generate a sequence of Xs
  x <- seq(max(0, floor(min(SMR.per.site[, Y]/1.1))), ceiling(max(SMR.per.site[, Y]) * 1.1), by=1)
  
  # Calculate the density of the function
  y.norm <- dnorm(x, mean=fit.norm$estimate[1],  sd=fit.norm$estimate[1])
  y.pois <- dpois(x, lambda=fit.pois$estimate[1])
  y.nbinom <- dnbinom(x, size=fit.nbinom$estimate[1], mu=fit.nbinom$estimate[2])
  
  # Store the result in the data.frame
  densities <- rbind(densities, data.frame(Y = rep(Y, length(x)), x=x, y.norm=y.norm, y.pois=y.pois, y.nbinom=y.nbinom))
  
  # Store the AIC
  AIC <- rbind(AIC, data.frame(Y=Y, AIC.norm=fit.norm$aic, AIC.pois=fit.pois$aic, AIC.nbinom=fit.nbinom$aic))

# End of the loop for the Y var
}

# Remove the first lines in the storage dataset
densities <- densities[-1, ]
AIC <- AIC[-1, ]

# Melt the data.frame for plotting purpose
SMR.per.site.melt <- melt(as.data.frame(SMR.per.site), measure.vars=Y.all)

# Melt the densities data.frame for plotting purpose
densities.melt <- melt(densities, id.vars=c('Y', 'x'), measures.vars=c('y.norm', 'y.pois', 'y.nbinom'))
names(densities.melt) <- c('variable', 'x', 'density', 'y')

# Multiple plotting with ggplot
pdf('out/data visualisation/tmp/Y_distributions.pdf', width=10, height=5)
ggplot(data=SMR.per.site.melt, aes(x=value)) +
  geom_histogram(aes(group=variable, y=..density..,), bins = 10) +
  geom_line(data=densities.melt, aes(x=x, y=y, color=density)) +
  facet_wrap(~variable, nrow=1, ncol=2, scales = 'free') +
  ggtitle('Distributions of the possible Y variables')
dev.off()

# Calculte best distribution for each Y variable
sapply(1:length(Y.all), function(w) which.min(AIC[w, 2:4])) #nbinomial for both

#### Part 4 : Generate functional observation of the explanatory variable per site

# Count the number of sites
nSites <- length(unique(SMR.per.site$Site))

# Create a data.frame to store the funobs
SMR.funobs <- as.data.frame(matrix(ncol=nSites + 3, nrow = length(Vars) * nKnots))
names(SMR.funobs) <- c('Var', 'Var.x', 'FunobsMean', paste0('Funobs', 1:nSites))

# Create the variable name
SMR.funobs$Var <- rep(Vars, each=nKnots)

# Find the range of the variable
SMR.Depth.range <- findRange(SMR.clean$Depth, 0)
SMR.Velocity.range <- findRange(SMR.clean$Velocity, 0)
SMR.D50.range <- findRange(SMR.clean$D50, 1)
SMR.D84.range <- findRange(SMR.clean$D84, 1)
SMR.Temp.range <- findRange(SMR.clean$Temp)

# Calculate the x values for this range
SMR.funobs$Var.x[which(SMR.funobs$Var=='Depth')] <- seq(SMR.Depth.range[1], SMR.Depth.range[2], length.out = nKnots)
SMR.funobs$Var.x[which(SMR.funobs$Var=='Velocity')] <- seq(SMR.Velocity.range[1], SMR.Velocity.range[2], length.out = nKnots)
SMR.funobs$Var.x[which(SMR.funobs$Var=='D50')] <- seq(SMR.D50.range[1], SMR.D50.range[2], length.out = nKnots)
SMR.funobs$Var.x[which(SMR.funobs$Var=='D84')] <- seq(SMR.D84.range[1], SMR.D84.range[2], length.out = nKnots)
SMR.funobs$Var.x[which(SMR.funobs$Var=='Temp')] <- seq(SMR.Temp.range[1], SMR.Temp.range[2], length.out = nKnots)

# Calculate the mean functional observation
SMR.funobs$FunobsMean[which(SMR.funobs$Var=='Depth')] <- density(SMR.clean$Depth, from = SMR.Depth.range[1], to = SMR.Depth.range[2], n = nKnots)$y
SMR.funobs$FunobsMean[which(SMR.funobs$Var=='Velocity')] <- density(SMR.clean$Velocity, from = SMR.Velocity.range[1], to = SMR.Velocity.range[2], n = nKnots)$y
SMR.funobs$FunobsMean[which(SMR.funobs$Var=='D50')] <- density(SMR.clean$D50, from = SMR.D50.range[1], to = SMR.D50.range[2], n = nKnots)$y
SMR.funobs$FunobsMean[which(SMR.funobs$Var=='D84')] <- density(SMR.clean$D84, from = SMR.D84.range[1], to = SMR.D84.range[2], n = nKnots)$y
SMR.funobs$FunobsMean[which(SMR.funobs$Var=='Temp')] <- density(SMR.clean$Temp, from = SMR.Temp.range[1], to = SMR.Temp.range[2], n = nKnots)$y

# Calculate the functional observation at each site
for (Site.i in 1:nSites) {
  SMR.funobs[which(SMR.funobs$Var=='Depth'), paste0('Funobs', Site.i)] <- density(SMR.clean$Depth[which(SMR.clean$Site == Site.i)], from = SMR.Depth.range[1], to = SMR.Depth.range[2], n = nKnots)$y
  SMR.funobs[which(SMR.funobs$Var=='Velocity'), paste0('Funobs', Site.i)] <- density(SMR.clean$Velocity[which(SMR.clean$Site == Site.i)], from = SMR.Velocity.range[1], to = SMR.Velocity.range[2], n = nKnots)$y
  SMR.funobs[which(SMR.funobs$Var=='D50'), paste0('Funobs', Site.i)] <- density(SMR.clean$D50[which(SMR.clean$Site == Site.i)], from = SMR.D50.range[1], to = SMR.D50.range[2], n = nKnots)$y
  SMR.funobs[which(SMR.funobs$Var=='D84'), paste0('Funobs', Site.i)] <- density(SMR.clean$D84[which(SMR.clean$Site == Site.i)], from = SMR.D84.range[1], to = SMR.D84.range[2], n = nKnots)$y
  SMR.funobs[which(SMR.funobs$Var=='Temp'), paste0('Funobs', Site.i)] <- density(SMR.clean$Temp[which(SMR.clean$Site == Site.i)], from = SMR.Temp.range[1], to = SMR.Temp.range[2], n = nKnots)$y
}

# Save the functional observations
saveRDS(SMR.funobs, 'data/SMR_2017_clean_fun_obs.Rda')

# Some visualisation, but first melting the data.frame
SMR.funobs.melt <- melt(SMR.funobs, id.vars = c('Var', 'Var.x'))

# Plotting with ggplot
pdf('out/data visualisation/tmp/SMR_funobs_all.pdf')
for (var in Vars) {
print(ggplot(data=SMR.funobs.melt[which(SMR.funobs.melt$Var==var), ], aes(x=Var.x, y=value)) +
      geom_line(aes(color=variable), show.legend=FALSE) +
      facet_wrap(~variable, nrow=5, ncol=6) +
      ylab(var) +
      ggtitle('Distribution at each site SMR', var))
}
dev.off()

#### Part 5 : We will focus on the 4 first site for visualisation of the function observation

# The vars we want to focus on 
var.list <- c('Depth', 'Velocity', 'D50', 'Temp')

# Function to transform variable to full names
TransVarName <- function(var)
  sapply(var, function(w) c('Depth (cm)', 'Velocity (m/s)', 'D50 (mm)', 'Temperature (ºC)')[var.list == w])

# Function to transform variable to full names
TransVarName2 <- function(var)
  sapply(var, function(w) c('a) Depth', 'b) Velocity', 'c) Median substrate size (D50)', 'd) Temperature')[var.list == w])

# Create a data.frame with all the values for the histogram
Hist.data <- SMR.clean[, c('Site', var.list)] %>% 
             melt(id.vars='Site', variable_name = 'Var') %>%
             filter(Site %in% 1:4) %>%
             mutate(Site.name = paste0('Site ', Site), Var.name = TransVarName(Var)) %>%
             dplyr::rename(x = value)
              
# Create a data.frame with all the values for the densities
Density.data <-  SMR.funobs.melt %>% 
                 filter(Var %in% var.list,
                 variable %in% c('Funobs1', 'Funobs2', 'Funobs3', 'Funobs4')) %>%
                 mutate(Site = as.integer(substr(variable, 7, 8))) %>%
                 mutate(Site.name = paste0('Site ', Site), Var.name = TransVarName(Var)) %>%
                 dplyr::rename(x = Var.x, y = value) %>%
                 dplyr::select(Site, Var, x, y, Site.name, Var.name)

# Create a data.frame with all the mean values per site
Mean.data <- Hist.data %>%
             group_by(Site, Var, Site.name, Var.name) %>%
             summarize(x = mean(x))

# We have to do some cleaning in the density.data
for (var.i in var.list) {
  for (site.i in 1:4) {
    factor <- ifelse (var.i == 'Temp', 1.05, 1.3)
    Limit <- c(max(min(Hist.data$x[Hist.data$Var == var.i & Hist.data$Site == site.i]) / factor, 0), max(Hist.data$x[Hist.data$Var == var.i & Hist.data$Site == site.i]) * factor)
    Density.data <- Density.data[-which(Density.data$Var == var.i & Density.data$Site == site.i & (Density.data$x<Limit[1] | Density.data$x>Limit[2])), ]
  }
}

# Generate the resulting plots for each variable
for (var.i in var.list)
{
  assign(paste0('p_', var.i), 
         ggplot(data=Hist.data[which(Hist.data$Var==var.i),], aes(x=x)) + 
           geom_histogram(aes(y=..density..),  bins=10, fill='#c1c1c1', color='#636363', lwd=0.25) +
           geom_line(data=Density.data[which(Density.data$Var==var.i),], aes(x=x, y=y, color='Kernel density estimate (KDE)       '), lwd=0.75) +
           geom_line(data=Density.data[which(Density.data$Var==var.i),], aes(x=x, y=y, color='Mean value'),  col=rgb(0,0,0,0)) +
           geom_vline(data=Mean.data[which(Mean.data$Var==var.i), ], aes(xintercept=x , color='Mean value'), lwd=0.75, linetype='dashed', show.legend=FALSE) +
           xlab(TransVarName(var.i)) + 
           ylab('') +
           labs(color='') + 
           scale_color_manual(values=c('#185fb1', '#ca3a27')) +
           facet_wrap(~Site.name, ncol=4, scales='free') +
           theme(axis.text=element_text(size=8),
                 axis.title=element_text(size=8))
  )
}

# Create a big figure with all the variables
fig <- ggarrange(p_Depth, p_Velocity, p_D50, p_Temp, 
                 nrow=4, 
                 labels = TransVarName2(var.list),
                 font.label = list(size=10, face='plain'),
                 vjust=-0.2,
                 hjust=0,
                 common.legend = TRUE,
                 legend='bottom')

# Generate the final .pdf version
pdf('out/data visualisation/Figure_3_Fun_obs.pdf')
annotate_figure(fig,
                top = text_grob(''),
                left = text_grob("Probability Distribution Function (PDF)", rot = 90, size=9, hjust = )
                )
dev.off()
