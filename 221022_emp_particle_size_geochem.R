################################################################################
#October 22-2022 Loon Lake Data Analysis
################################################################################
#Loon Lake Empirical Mode Decomposition 
################################################################################

###################################################
#LOADING LIBRARIES AND SETTING WORKING DIRECTORIES
###################################################

#Loading libraries

#Run for the first time only:
# install.packages("librarian")
# librarian::lib_startup(librarian, global = TRUE)
# source: https://cran.r-project.org/web/packages/librarian/readme/README.html

librarian::shelf(ggplot2,dplyr,imputeTS,curl,EMD,plot3D,plot3Drgl)


# library(ggplot2)
# library(reshape2)
# library(lubridate)
# library(gridExtra)
# library(grid)
# library(plyr)
# library(dplyr)
# library(nlme)
# library(doBy)
# library(MASS)
# library(lsmeans)
# library(carData)
# library(utils)
# library(multcompView)
# library(cowplot)
# library(scales)
# library(imputeTS)
# library(RColorBrewer)
# library(tidyr)
# library(gridExtra)
# library(lattice)
# library(grid)
# library(gridBase)
# library(MuMIn)
# library(hexbin)
# library(MASS)
# library(devtools)
# library(curl)
# library(scales)
# library(ggpubr)
# library(psych)
# library(tidyr)
# library(EMD)



#Using Hatten's theme for plots

theme_httn<-  theme(axis.text=element_text(colour="black",size=14),
                    axis.title = element_text(size = 20),
                    panel.grid.minor= element_line(colour = "gray", linetype = "dotted"), 
                    panel.grid.major = element_line(colour = "gray", linetype = "dashed"),
                    panel.border = element_rect(fill=NA, colour = "black", size=1),
                    panel.background=element_rect(fill="white"),
                    axis.ticks.length = unit(0.254, "cm"),
                    axis.ticks = element_line(colour = "black", size=1), 
                    axis.line = element_line(colour = "black"),
                    legend.position = c(0.92,0.85),
                    legend.direction = "vertical",
                    legend.background = element_blank(),
                    legend.key.size = unit(1.0, 'lines'),#Changing spacing between legend keys
                    legend.title = element_text())

#creating a theme for a basic plot
theme_bp<-theme(axis.title.x = element_text(size=14),
                axis.title.y = element_text(size=14),
                axis.line.x = element_line(size=1),
                axis.ticks.x = element_line(size=1),
                axis.line.y = element_line(size=1),
                axis.ticks.y = element_line(size=1),
                axis.text.x = element_text(size=12),
                axis.text.y = element_text(size=12))


###################################################
#Loading datasets
##################################################

eda<-read.csv("180912_eda_gls_results_lint.csv")
eda$er_n<-factor((eda$er_n),levels = rev(c("No-Event","1996 A.D.","1982 A.D.","1971 A.D.","1964-65 A.D.","1946 A.D.","1890 A.D.","1770 A.D.","1690 A.D.","1600 A.D.","1470 A.D.",
"1370 A.D.","1350 A.D.","1260 A.D.","1240 A.D.","1150 A.D.","1030 A.D.","940 A.D.","920 A.D.","760 A.D.","740 A.D.","710 A.D.","690 A.D.","660 A.D.","630 A.D.", "550 A.D.","510 A.D.")))

eda$prd_c<-factor((eda$prd_c),levels = (c("400-500","500-600","600-700","700-800","800-900","900-1000",
                                          "1000-1100","1100-1200","1200-1300","1300-1400","1400-1500",
                                          "1500-1600","1600-1700","1700-1800","1800-1900","1900-present")))

#################################################
#Entropy

### Empirical Mode Decomposition (hrel_6)
par(mfrow=c(3,1), mar=c(1.5,0.1,0.2,0.1))
im_hrel <- emd(eda$hrel_6, eda$z, boundary="wave")


### Ploting the IMF's
par(mfrow=c(3,1), mar=c(1.5,0.1,0.2,0.1))
par(mfrow=c(im_hrel$nimf+1, 1), mar=c(1.5,0.1,0.2,0.1))
rangeimf <- range(im_hrel$imf)
for(i in 1:im_hrel$nimf)
   plot(eda$z, im_hrel$imf[,i], type="l", xlab="",
          ylab="", ylim=rangeimf, main=
            paste(i, "-th IMF", sep="")); abline(h=0)
plot(eda$z, im_hrel$residue, xlab="", ylab="",
        main="residue", type="l")


### Empirical Mode Decomposition (d_50)
par(mfrow=c(3,1), mar=c(0.5,0.1,0.2,0.1))
try <- emd(eda$d50, eda$z, boundary="wave")


### Ploting the IMF's
par(mfrow=c(3,1), mar=c(0.5,0.1,0.2,0.1))
par(mfrow=c(try$nimf+1, 1), mar=c(0.5,0.1,0.2,0.1))
rangeimf <- range(try$imf)
for(i in 1:try$nimf)
  plot(eda$z, try$imf[,i], type="l", xlab="",
       ylab="", ylim=rangeimf, main=
         paste(i, "-th IMF", sep="")); abline(h=0)
plot(eda$z, try$residue, xlab="", ylab="",
     main="residue", type="l")

#Apparently the event c.a. 690 A.D. (which we suspect is from local nature)
#eclipses most of the other changes across the particle size distribution.
#Let's try expressing those changes in d50, using the phi-scale ~log10

par(mfrow=c(3,1), mar=c(0.5,0.1,0.2,0.1))
try <- emd(log10(eda$d50), eda$z, boundary="wave")


### Ploting the IMF's
par(mfrow=c(3,1), mar=c(0.5,0.1,0.2,0.1))
par(mfrow=c(try$nimf+1, 1), mar=c(0.5,0.1,0.2,0.1))
rangeimf <- range(try$imf)
for(i in 1:try$nimf)
  plot(eda$z, try$imf[,i], type="l", xlab="",
       ylab="", ylim=rangeimf, main=
         paste(i, "-th IMF", sep="")); abline(h=0)
plot(eda$z, try$residue, xlab="", ylab="",
     main="residue", type="l")

#Silt
par(mfrow=c(3,1), mar=c(0.5,0.1,0.2,0.1))
try <- emd(log10(eda$silt), eda$z, boundary="wave")


### Ploting the IMF's
par(mfrow=c(3,1), mar=c(0.5,0.1,0.2,0.1))
par(mfrow=c(try$nimf+1, 1), mar=c(0.5,0.1,0.2,0.1))
rangeimf <- range(try$imf)
for(i in 1:try$nimf)
  plot(eda$z, try$imf[,i], type="l", xlab="",
       ylab="", ylim=rangeimf, main=
         paste(i, "-th IMF", sep="")); abline(h=0)
plot(eda$z, try$residue, xlab="", ylab="",
     main="residue", type="l")

#############
#Geochemical Indicators

#oc
par(mfrow=c(3,1), mar=c(0.5,0.1,0.2,0.1))
try <- emd((eda$oc), eda$z, boundary="wave")


### Ploting the IMF's
par(mfrow=c(3,1), mar=c(0.5,0.1,0.2,0.1))
par(mfrow=c(try$nimf+1, 1), mar=c(0.5,0.1,0.2,0.1))
rangeimf <- range(try$imf)
for(i in 1:try$nimf)
  plot(eda$z, try$imf[,i], type="l", xlab="",
       ylab="", ylim=rangeimf, main=
         paste(i, "-th IMF", sep="")); abline(h=0)
plot(eda$z, try$residue, xlab="", ylab="",
     main="residue", type="l")

#d13c
par(mfrow=c(3,1), mar=c(0.5,0.1,0.2,0.1))
try <- emd((eda$d13ci), eda$z, boundary="wave")


### Ploting the IMF's
par(mfrow=c(3,1), mar=c(0.5,0.1,0.2,0.1))
par(mfrow=c(try$nimf+1, 1), mar=c(0.5,0.1,0.2,0.1))
rangeimf <- range(try$imf)
for(i in 1:try$nimf)
  plot(eda$z, try$imf[,i], type="l", xlab="",
       ylab="", ylim=rangeimf, main=
         paste(i, "-th IMF", sep="")); abline(h=0)
plot(eda$z, try$residue, xlab="", ylab="",
     main="residue", type="l")

#d15n
par(mfrow=c(3,1), mar=c(0.5,0.1,0.2,0.1))
try <- emd((eda$d15ni), eda$z, boundary="wave")


### Ploting the IMF's
par(mfrow=c(3,1), mar=c(0.5,0.1,0.2,0.1))
par(mfrow=c(try$nimf+1, 1), mar=c(0.5,0.1,0.2,0.1))
rangeimf <- range(try$imf)
for(i in 1:try$nimf)
  plot(eda$z, try$imf[,i], type="l", xlab="",
       ylab="", ylim=rangeimf, main=
         paste(i, "-th IMF", sep="")); abline(h=0)
plot(eda$z, try$residue, xlab="", ylab="",
     main="residue", type="l")
