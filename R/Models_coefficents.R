####################################################################################
#     Modelling fish habitat selection using functional regression 🐟
#     Created by : Jeremie Boudreault (Jeremie.Boudreault@ete.inrs.ca)
#     Current file : Models_Coefficients.R
#     Aim : Extract and plot the models coefficients for GAMs and FRMs
####################################################################################

#### Part 0 : Libraries, functions and variable definition

# Import library
library('dplyr') # Data.frame manipulation
library('ggplot2') # Data visualisation
library('mgcv') # For GAMs
library('FDboost') # For FRMs
library('ggpubr')
library('RColorBrewer') # For color of the plot

# Import SMR per site values
SMR.per.site <- readRDS('data/SMR_2017_clean_per_site.Rda')

# Import GAMs models
GAM.fry <- readRDS('out/models/GAM_nFry_M.rda')
GAM.parr1 <- readRDS('out/models/GAM_nParr1_M.rda')

# Import FRMs models
FRM.fry <- readRDS('out/models/FRM_nFry_M.rda')
FRM.parr1 <- readRDS('out/models/FRM_nParr1_M.rda')

# Create the color for the plots
myColors <- brewer.pal(9, "Set1")[c(2, 3, 4, 1)]
names(myColors) <- c('Depth', 'Velocity', 'D50', 'Temperature')
colScale <- scale_colour_manual(name = "grp",values = myColors)
fillScale <- scale_fill_manual(name = "grp",values = myColors)

##### Part 1: Visualisation for FRMs

# First extract values of the coefficient for the fry 
coef(FRM.fry) # Order is velocity, D50 and Temp
CoefTmp <- sapply(1:3, function(w) coef(FRM.fry)$smterms[[w]]$value) # Extract the coefs
CoefTmp.X <- sapply(1:3, function(w) coef(FRM.fry)$smterms[[w]]$x) # Extract the x values
FRM.fry.coef <- data.frame(var=rep(c('Velocity', 'D50', 'Temp'), each=nrow(CoefTmp.X)),
                           x = as.vector(CoefTmp.X),
                           y = as.vector(CoefTmp))

# Import the leave-one-out coefficient for estimation of confidence interval (Fry)
FRM.fry.coef.loo <- readRDS('out/coefficients/Coef_LOO_FRM_nFry_M.rda')
FRM.fry.coef.loo <- FRM.fry.coef.loo[which(FRM.fry.coef.loo$var %in% c('Velocity', 'D50', 'Temp')), ]
FRM.fry.coef.loo$sd <- apply(FRM.fry.coef.loo[, 3:28], MARGIN = 1, FUN=sd)
FRM.fry.coef.loo <- FRM.fry.coef.loo[, -(3:28)]

# Combine the two tables (Estimates + Leave-one-out)
FRM.fry.coef <- FRM.fry.coef %>% left_join(FRM.fry.coef.loo, by=c('x', 'var'))

# Generate some NA for plotting purposes
FRM.fry.coef <- rbind(FRM.fry.coef, data.frame(var='Depth', x = c(0, 30, 103, NA, NA, NA), y=c(NA, NA, NA, -2, 0, 1.4), sd=rep(0, 6)))
FRM.fry.coef$var <- factor(FRM.fry.coef$var, levels = c('Depth', 'Velocity', 'D50', 'Temperature'))
FRM.fry.coef$var[is.na(FRM.fry.coef$var)] <- 'Temperature'

# Plot the resulting coefficients for Fry
p_fry <- ggplot(data=FRM.fry.coef, aes(x=x, y=y)) +
  geom_ribbon(aes(ymin=y-1.96*sd, ymax=y+1.96*sd, x=x, fill=var), alpha=0.15,show.legend=FALSE)+
  geom_line(aes(color=var), lwd=0.7, show.legend=FALSE) +
  geom_hline(aes(yintercept=0), color='black', lty=2, lwd=0.5, alpha=0.3) + 
  facet_wrap(~var, nrow=1, ncol=4, scales='free') +
  colScale + fillScale +
  xlab('') + ylab('') +
  theme(axis.title.y = element_blank(),
        axis.text=element_text(size=8))
p_fry

# Second, extract values of the coefficient for the 1+ parr
coef(FRM.parr1) # Order is Depth, D50 and Temp
CoefTmp <- sapply(1:3, function(w) coef(FRM.parr1)$smterms[[w]]$value) # Extract the coefs
CoefTmp.X <- sapply(1:3, function(w) coef(FRM.parr1)$smterms[[w]]$x) # Extract the x values
FRM.parr1.coef <- data.frame(var=rep(c('Depth', 'D50', 'Temp'), each=nrow(CoefTmp.X)),
                           x = as.vector(CoefTmp.X),
                           y = as.vector(CoefTmp))

# Import the leave-one-out coefficient for estimation of confidence interval (1+ Parr)
FRM.parr1.coef.loo <- readRDS('out/coefficients/Coef_LOO_FRM_nParr1_M.rda')
FRM.parr1.coef.loo <- FRM.parr1.coef.loo[which(FRM.parr1.coef.loo$var %in% c('Depth',  'D50', 'Temp')), ]
FRM.parr1.coef.loo$sd <- apply(FRM.parr1.coef.loo[, 3:28], MARGIN = 1, FUN=sd)
FRM.parr1.coef.loo <- FRM.parr1.coef.loo[, -(3:28)]

# Combine the two tables (Estimates + Leave-one-out)
FRM.parr1.coef <- FRM.parr1.coef %>% left_join(FRM.parr1.coef.loo, by=c('x', 'var'))

# Generate some NA for plotting purposes
FRM.parr1.coef <- rbind(FRM.parr1.coef, data.frame(var=rep('Velocity', 6), x = c(0, 1, 1.9, NA, NA, NA), y=c(NA, NA, NA, 1.25, 0, -9.5), sd=rep(0, 6)))
FRM.parr1.coef$var <- factor(FRM.parr1.coef$var, levels = c('Depth', 'Velocity', 'D50', 'Temperature'))
FRM.parr1.coef$var[is.na(FRM.parr1.coef$var)] <- 'Temperature'

# Plot the resulting coefficients for 1+ parr
p_parr1 <- ggplot(data=FRM.parr1.coef, aes(x=x, y=y)) +
  geom_ribbon(aes(ymin=y-1.96*sd, ymax=y+1.96*sd, x=x, fill=var), alpha=0.15, show.legend=FALSE)+
  geom_line(aes(color=var), lwd=0.7, show.legend=FALSE) +
  geom_hline(aes(yintercept=0), color='black', lty=2, lwd=0.5, alpha=0.3) + 
  facet_wrap(~var, nrow=1, ncol=4, scales='free') +
  colScale + fillScale +
  xlab('') + ylab('') +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text=element_text(size=8))
p_parr1

# Generate a common figure
fig <- ggarrange(p_fry, p_parr1, 
                 nrow=2,
                 labels = c('a) Fry (0+)', 'b) Parr (1+)'),
                 font.label = list(size=10, face='plain'),
                 vjust=-0.4,
                 hjust=-0.2,
                 common.legend = TRUE,
                 legend='bottom')

# Generate the final figure for FRM
pdf('out/coefficients/Figure_4_FRMs_Coefficients.pdf', height = 5, width=8)
annotate_figure(fig,
                top = text_grob(''),
                left = text_grob("Models coefficients", rot = 90, size=10),
                bottom = text_grob("Depth (cm)                                           Velocity (m/s)                                    D50 (mm)                               Temperature (ºC)", size=9, hjust=0.48, vjust=-0.3))
dev.off()

###### Part 2 : Visualisation for the GAMs

# First, extract the coefficient for Fry
x.temp <- data.frame(Temp = seq(floor(min(SMR.per.site$Temp)), ceiling(max(SMR.per.site$Temp)), length.out=40))
coef <- predict(GAM.fry,type = 'terms', se.fit=T, newdata = data.frame(x.temp))
GAM.fry.coef <- data.frame(var='Temp', 
                           x=x.temp[, 1],
                           y=as.vector(coef$fit[, 1]),
                           se=as.vector(coef$se.fit[, 1]))

# Generate NA values and arrange factors for plotting purpopose
GAM.fry.coef <- rbind(GAM.fry.coef, data.frame(var=rep(c('Depth', 'D50'), each=4),  x=c(25, 70, NA, NA, 75, 175, NA, NA), y=c(NA, NA, -0.5, 1, NA, NA, -1, 0.5), se=rep(0, 8)))
GAM.fry.coef$var <- factor(GAM.fry.coef$var, levels = c('Depth', 'Velocity', 'D50', 'Temperature'))
GAM.fry.coef$var[is.na(GAM.fry.coef$var)] <- 'Temperature'

# Generate the plot for Fry
p_fry <- ggplot(data=GAM.fry.coef, aes(x=x, y=y)) +
  geom_ribbon(aes(ymin=y-se, ymax=y+se, x=x, fill=var), alpha=0.15,show.legend=FALSE)+
  geom_line(aes(color=var), lwd=0.7, show.legend=FALSE) +
  geom_hline(aes(yintercept=0), color='black', lty=2, lwd=0.5, alpha=0.3) + 
  facet_wrap(~var, nrow=1, ncol=4, scales='free') +
  colScale + fillScale +
  xlab('') + ylab('') +
  theme(axis.title.y = element_blank(),
        axis.text=element_text(size=8))
p_fry

# Second, extract the coefficient for Parr1
x.temp <- data.frame(Depth = seq(min(SMR.per.site$Depth), max(SMR.per.site$Depth), length.out=40),
                     D50 = seq(min(SMR.per.site$D50), max(SMR.per.site$D50), length.out=40),
                     Temp = seq(floor(min(SMR.per.site$Temp)), ceiling(max(SMR.per.site$Temp)), length.out=40))
coef <- predict(GAM.parr1, type = 'terms', se.fit=T, newdata = x.temp)
GAM.parr1.coef <- data.frame(var=rep(c('Depth', 'D50', 'Temperature'), each=40),
                             x=as.vector(unlist(x.temp)),
                             y=as.vector(unlist(coef$fit)),
                             se=as.vector(unlist(coef$se.fit)))
GAM.parr1.coef$var <- factor(GAM.parr1.coef$var, levels = c('Depth', 'Velocity', 'D50', 'Temperature'))


# Plot the resulting coefficients for 1+ parr
p_parr1 <- ggplot(data=GAM.parr1.coef, aes(x=x, y=y)) +
  geom_ribbon(aes(ymin=y-se, ymax=y+se, x=x, fill=var), alpha=0.15, show.legend=FALSE)+
  geom_line(aes(color=var), lwd=0.7, show.legend=FALSE) +
  geom_hline(aes(yintercept=0), color='black', lty=2, lwd=0.5, alpha=0.3) + 
  facet_wrap(~var, nrow=1, ncol=4, scales='free') +
  colScale + fillScale +
  xlab('') + ylab('') +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text=element_text(size=8))
p_parr1

# Generate a common figure
fig <- ggarrange(p_fry, p_parr1, 
                 nrow=2,
                 labels = c('a) Fry (0+)', 'b) Parr (1+)'),
                 font.label = list(size=10, face='plain'),
                 vjust=-0.4,
                 hjust=-0.2,
                 common.legend = TRUE,
                 legend='bottom')

# Generate the final figure for GAM
pdf('out/coefficients/Figure_5_GAMs_Coefficients.pdf', height = 5, width=7)
annotate_figure(fig,
                top = text_grob(''),
                left = text_grob("Models coefficients", rot = 90, size=10),
                bottom = text_grob("Depth (cm)                                                  D50 (mm)                                         Temperature (ºC)", size=9, hjust=0.45, vjust=-0.3))
dev.off()

###### Part 3 : Comparison between GAMs and FRMs for D50 and 1+ parr

# Extract the D50 coefficients for 1+ parr (FRM)
FRM.D50 <- FRM.parr1.coef %>%
           dplyr::filter(var == 'D50')  %>%                # Select D50
           dplyr::mutate(se = sd * 1.96, model='FRM')  %>% # Change standard deviation to std error (including the 95% CI)
           dplyr::select(-sd)                           # Remove STD
  
# Extract the D50 coefficients for 1+ parr (GAM)
GAM.D50 <- GAM.parr1.coef %>%
           dplyr::filter(var == 'D50') %>%
           dplyr::mutate(model='GAM')

# Extract the extrapolate values for the GAM
x.temp <- data.frame(Depth = seq(min(SMR.per.site$Depth), max(SMR.per.site$Depth), length.out=40),
                     D50 = FRM.D50$x,
                     Temp = seq(floor(min(SMR.per.site$Temp)), ceiling(max(SMR.per.site$Temp)), length.out=40))
coef <- predict(GAM.parr1, type = 'terms', se.fit=T, newdata = x.temp)

# Generate a data frame for GAM (Extrapolate)
GAM.D50.ex <- data.frame(var=rep('D50', 40), x=FRM.D50$x, y=coef$fit[, 2], se=coef$se.fit[, 2], model='GAM (extrapolate)') %>%
                       mutate(ex.lower = (x <= GAM.D50$x[3]), 
                              ex.upper = (x >= GAM.D50$x[length(GAM.D50$x)]))

myColors <- brewer.pal(9, "Set1")[c(4, 4, 1)]
colScale <- scale_colour_manual(values = myColors)

# Plot the first graph (FRM)
p_frm <- 
ggplot(data=FRM.D50, aes(x=x, y=y)) +
  geom_line(aes(linetype='Extrapolated', color='Extrapolated')) +
  geom_line(aes(linetype='Extracted coefficient', color='Extracted coefficient')) +
  geom_line(color=myColors[1], lwd=1.1) +
  geom_vline(data=data.frame(x=c(min(GAM.D50$x), max(GAM.D50$x))), aes(xintercept=x), alpha=0.4, lwd=0.5, lty=3) +
  geom_ribbon(aes(ymin=y-se, ymax=y+se), fill=myColors[1], alpha=0.15)  +
  scale_colour_manual("", values=myColors) +
  scale_linetype_manual("", values=c('solid', 'dashed')) +
  xlab('') + ylab('') +
  scale_y_continuous(limits=c(-2.6, 0.7)) +
  scale_x_continuous(limits=c(0, 310)) +
  theme(axis.title.y = element_blank(),
        axis.text=element_text(size=8))

# Plot the second graph (GAM)
p_gam <- 
ggplot(data=GAM.D50, aes(x=x, y=y)) +
  geom_line(data=GAM.D50.ex[GAM.D50.ex$ex.upper, ], aes(x=x, y=y), color=myColors[1], lwd=1.1, linetype=2) +
  geom_line(data=GAM.D50.ex[GAM.D50.ex$ex.lower, ], aes(x=x, y=y), color=myColors[1], lwd=1.1, linetype=2) +
  geom_line(color=myColors[1], lwd=1.1) +
  geom_vline(data=data.frame(x=c(min(GAM.D50$x), max(GAM.D50$x))), aes(xintercept=x), alpha=0.4, lwd=0.5, lty=3) +
  geom_ribbon(data=GAM.D50.ex, aes(x=x, ymin=y-se, ymax=y+se), fill=myColors[1], alpha=0.15) +
  xlab('D50 (mm)') + ylab('') +
  scale_y_continuous(limits=c(-2.6, 0.7)) +
  scale_x_continuous(limits=c(0, 310)) +
  theme(axis.title.y = element_blank(),
        axis.title.x=element_text(size=9),
        axis.text=element_text(size=8))

# Generate a common figure
fig <- ggarrange(p_frm, p_gam, 
                 nrow=2,
                 labels = c('a) Functional regression model (FRM)', 'b) Generalized additive model (GAM)'),
                 font.label = list(size=10, face='plain'),
                 vjust=-0.4,
                 hjust=-0.15,
                 common.legend = TRUE,
                 legend='bottom')


# Generate the final figure for GAM
pdf('out/coefficients/Figure_6_D50_Parr1_Comparison.pdf', height = 6, width=4)
annotate_figure(fig,
                top = text_grob(''),
                left = text_grob("Models coefficients", rot = 90, size=10))
dev.off()


