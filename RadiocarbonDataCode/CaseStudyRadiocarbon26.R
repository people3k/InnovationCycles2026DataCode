
#Script for analyzing the radiocarbon data for the continental US and American Southwest archaeological regions associated with
#Freeman et al. 2026: Long-term cycles in food producing archaeological regions

##Load Packages==================================
library(elevatr)
library(terra)
library(ggplot2)
library(dplyr)
library(zoo)
library(rcarbon)
library(cowplot)
library(tidyverse)
library(sf)
library(maps)
library(rnaturalearth)
library(rnaturalearthdata)

##Set working directory

##1. Run KDEs for each case study
####Sonoran Desert, Case ID 5 ===========================================================
SPD<-read.csv(file="data/RawP3Kc14.csv", header=T)
boxsd<- subset(SPD, Latitude>23 & Latitude<35 & Longitude>-115 & Longitude< -108)
#write.table(boxsd, file = "NERDv4_0/SonoranDesert.csv", sep = ",", col.names=NA)
#boxsd<-read.csv(file="data/SonoranDesert.csv", header=T)
###MAP
counties<-map_data("state")

ArchGlobeMap<-ggplot() +
  geom_polygon(data = counties, mapping = aes(x = long, y = lat, group = group),
               fill = "grey", color = "white") +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45)+
  #geom_polygon(data = canada, aes(x=long, y = lat, group = group),
  #    fill = "white", color="black") +
  theme_bw()+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Longitude", y="Latitude", title = "US ArchaeoGlobe Regions and Radiocarbon")+
  geom_point(data=boxsd, aes(Longitude, Latitude, color=factor(Province)),
             inherit.aes = FALSE, alpha = 0.5, size = 2)
ArchGlobeMap
#remove NA's from the siteID column
boxsd <- boxsd[!is.na(boxsd$SiteID), ]

CalMz <- calibrate(x = boxsd$Age,  errors = boxsd$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = boxsd$SiteID, ages = boxsd$Age, h = 100)

####Run SPD
spd.mz <- spd(CalMz, bins=boxbins, runm=200, timeRange=c(8200,200))
plot(spd.mz, runm=200, xlim=c(8200,200), type="simple")

calBP<-spd.mz$grid$calBP
PrDens<-spd.mz$grid$PrDens

##Check the effect of h function on SPD if you desire
#binsense(x=CalMz,y=SPD$SiteName,h=seq(0,500,100),timeRange=c(4000,100))

##KDE
####make KDEs
US.randates = sampleDates(CalMz, bins=boxbins, nsim=500,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(8200,200),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')
D.ckde$timeRange

##Write matrix of KDEs as a data frame
Check<-as.data.frame(D.ckde$res.matrix)
#I then convert NAs at the end due to KDE bandwidth to zeros to easily cbind the data frames below
Check2 <- replace(Check, is.na(Check), 0)
#Calculate the mean KDE of each time step of the 200 simulations
MKDE<-rowMeans(Check2)
#calculate the 5th and 95th percentile of each time step of the 200 KDE simulations
lo<-apply( Check2, #select columns
           1, # row-wise calcs
           quantile, probs=0.05) # give `quantile`
hi <-apply( Check2, #select columns
            1, # row-wise calcs
            quantile, probs=0.95) # give `quantile`

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<-dd %>%  filter(MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/SonoranKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/SonoranKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........
library(zoo)
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')
out200<-rollapply(dd2c,200,(sum),by=200,by.column=TRUE,align='right')
###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)

calBP<-c(8170, 8140, 8110, 8080, 8050, 8020, 7990, 7960, 7930, 7900,7870,7840,7810,7780,7750,7720,
         7690,7660,7630,7600,7570,7540,7510,7480,7450,7420,7390,7360,7330,7300,7270,7240,7210,7180, 7150, 7120, 7090,
         7060, 7030, 7000, 6970, 6940, 6910, 6880, 6850, 6820, 6790, 6760, 6730, 6700, 6670, 6640,
         6610, 6580, 6550, 6520, 6490, 6460, 6430, 6400, 6370, 6340, 6310, 6280, 6250, 6220, 6190, 6160, 6130,6100,
         6070,6040,6010,5980,5950,5920,5890,5860,5830,5800,5770,5740,5710,5680,5650,5620,5590,5560,5530,
         5500,5470, 5440,5410, 5380,5350,5320,5290,5260, 5230,5200,5170,5140,5110,5080,5050,5020,4990, 4960,
         4930,4900,4870,4840,4810,4780,4750,4720,4690,4660,4630,4600,4570,4540,4510,4480, 4450,4420,4390,4360, 4330,4300,
         4270,4240,4210,4180, 4150,4120,4090,4060,4030,4000,3970,3940,3910,3880,3850,3820,3790,3760,3730,3700,
         3670,3640,3610,3580,3550,3520,3490,3460,3430,3400,3370, 3340, 3310,3280, 3250, 3220, 3190,
         3160, 3130, 3100, 3070, 3040, 3010, 2980, 2950,2920, 2890,2860,2830,2800, 2770,2740,2710, 2680,
         2650, 2620,2590,2560, 2530, 2500,2470,2440,2410, 2380,2350, 2320,2290, 2260, 2230, 2200, 2170, 2140,2110,2080,
         2050,2020,1990,1960,1930,1900,1870,1840,1810,1780,1750,1720,1690,1660,1630,1600,1570,1540,1510,
         1480,1450,1420,1390,1360,1330,1300,1270,1240,1210,1180,1150,1120,1090,1060,1030,1000,970,
         940,910,880,850,820,790,760,730,700,670,640,610,580,550,520,490,460,430,400,370,340,310,280,250, 220)
sums<-cbind(calBP, MKDE, out50)

write.table(sums, file = "data/Sumbin/SonoranSumbin.csv", sep = ",", col.names=NA)

##200 year bins
MKDE<-rowMeans(out200)
##Add in the 30 year bin dates
calBP<-c(8000, 7800, 7600, 7400, 7200, 7000, 6800, 6600, 6400, 6200, 6000, 5800, 5600, 5400,
         5200, 5000, 4800, 4600, 4400, 4200, 4000, 3800, 3600, 3400, 3200, 3000, 2800, 2600,
         2400, 2200, 2000, 1800, 1600, 1400, 1200, 1000, 800, 600, 400
)
sums<-cbind(calBP, MKDE, out200)
write.table(sums, file = "data/Sumbin/SonoranSumbin200.csv", sep = ",", col.names=NA)


###calculate growth rates for 200 year summed bins

d <-read_csv("data/Sumbin/SonoranSumbin200.csv") %>%
  dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))

pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200 ,700, 1780, 2800, 3130, 4050, 4500, 7500, 8200),
                         labels=c('Post-Hohokam','Hohokam', 'Cienega', 'San Pedro', 'Silver Bell','Late Archaic','Middle Archaic','Early Archaic'))

write.table(pcgrowth, file = "data/Percapita200/SonoranPerCap200.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps.

d <-read_csv("data/Sumbin/SonoranSumbin.csv") %>%
  dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))

##Add ID variable for culture history periods
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200 ,700, 1780, 2800, 4050, 7500, 8200),
                         labels=c('Post-Hohokam','Hohokam', 'Cienega', 'Early Agricultural Period','Middle Archaic','Early Archaic'))


write.table(pcgrowth, file = "data/Percapita/SonoranPerCap.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the Sonoran desert
son30pc<- read.csv("data/Percapita/SonoranPerCap.csv")

#son30pc<- read.csv("data/Percapita200/SonoranPerCap200.csv")

son30pc2<-subset(son30pc, calBP<4000 & calBP>200)

#Standardize the mean KDE by the maximum mean KDE during the Neolithic 
StKDE<-(son30pc2$MKDE-min(son30pc2$MKDE))/(max((son30pc2$MKDE)-min(son30pc2$MKDE)))
##Add the standardized KDE to the Neolithic dataframe
son30pc3<-cbind(StKDE, son30pc2)

##Write food production file
#write.table(son30pc3, file = "data/Percapita2/SonoranPerCap.csv", sep = ",", col.names=NA)

pcSon <- ggplot(son30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  scale_y_continuous(limits=c(-.3,0.3))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE density", y="KDE per capita growth", title = "B. Sonoran KDE Per Capita Growth vs. Density")
  #geom_vline(xintercept = 0.20)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcSon

Soncpt <- ggplot(son30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "A. Sonoran Density vs. Time")
#  geom_hline(yintercept = 0.20)
#geom_vline(xintercept = 2550)+
#geom_vline(xintercept = 2250)+
#geom_vline(xintercept = 830)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
Soncpt

Figsonora<-plot_grid(Soncpt, pcSon, ncol=1, align="hv", axis = "rl")
Figsonora

#pdf("data/Figs/Sonora.pdf", width=17.55, height=15)
#Figsonora
#dev.off()

#####################N. Colorado Plat. (Fremont) CASE ID 7=======================================================
SPD<-read.csv(file="data/RawP3Kc14.csv", header=T)
boxsd<- subset(SPD, Latitude>38 & Latitude<42 & Longitude>-115 & Longitude< -109)
#write.table(boxsd, file = "data/Fremontdates.csv", sep = ",", col.names=NA)
#boxsd<-read.csv(file="data/Fremontdates.csv", header=T)
###MAP
counties<-map_data("state")

ArchGlobeMap<-ggplot() +
  geom_polygon(data = counties, mapping = aes(x = long, y = lat, group = group),
               fill = "grey", color = "white") +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45)+
  #geom_polygon(data = canada, aes(x=long, y = lat, group = group),
  #    fill = "white", color="black") +
  theme_bw()+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Longitude", y="Latitude", title = "US ArchaeoGlobe Regions and Radiocarbon")+
  geom_point(data=boxsd, aes(Longitude, Latitude, color=factor(Province)),
             inherit.aes = FALSE, alpha = 0.5, size = 2)
ArchGlobeMap

boxsd <- boxsd[!is.na(boxsd$SiteID), ]

CalMz <- calibrate(x = boxsd$Age,  errors = boxsd$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = boxsd$SiteID, ages = boxsd$Age, h = 100)

####Run SPD
spd.mz <- spd(CalMz, bins=boxbins, runm=200, timeRange=c(8200,200))
plot(spd.mz, runm=200, xlim=c(8200,200), type="simple")

##Check the effect of h function on SPD if you desire
#binsense(x=CalMz,y=SPD$SiteName,h=seq(0,500,100),timeRange=c(4000,100))

##KDE
####make KDEs
US.randates = sampleDates(CalMz, bins=boxbins, nsim=500,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(8200,200),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')
#D.ckde$timeRange

##Write matrix of KDEs as a data frame
Check<-as.data.frame(D.ckde$res.matrix)
#I then convert NAs at the end due to KDE bandwidth to zeros to easily cbind the data frames below
Check2 <- replace(Check, is.na(Check), 0)
#Calculate the mean KDE of each time step of the 200 simulations
MKDE<-rowMeans(Check2)
#calculate the 5th and 95th percentile of each time step of the 200 KDE simulations
lo<-apply( Check2, #select columns
           1, # row-wise calcs
           quantile, probs=0.05) # give `quantile`
hi <-apply( Check2, #select columns
            1, # row-wise calcs
            quantile, probs=0.95) # give `quantile`

calBP<-spd.mz$grid$calBP
PrDens<-spd.mz$grid$PrDens

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<- subset(dd, MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/FremontKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/FremontKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')
out200<-rollapply(dd2c,200,(sum),by=200,by.column=TRUE,align='right')
###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
calBP<-c(8170, 8140, 8110, 8080, 8050, 8020, 7990, 7960, 7930, 7900,7870,7840,7810,7780,7750,7720,
         7690,7660,7630,7600,7570,7540,7510,7480,7450,7420,7390,7360,7330,7300,7270,7240,7210,7180, 7150, 7120, 7090,
         7060, 7030, 7000, 6970, 6940, 6910, 6880, 6850, 6820, 6790, 6760, 6730, 6700, 6670, 6640,
         6610, 6580, 6550, 6520, 6490, 6460, 6430, 6400, 6370, 6340, 6310, 6280, 6250, 6220, 6190, 6160, 6130,6100,
         6070,6040,6010,5980,5950,5920,5890,5860,5830,5800,5770,5740,5710,5680,5650,5620,5590,5560,5530,
         5500,5470, 5440,5410, 5380,5350,5320,5290,5260, 5230,5200,5170,5140,5110,5080,5050,5020,4990, 4960,
         4930,4900,4870,4840,4810,4780,4750,4720,4690,4660,4630,4600,4570,4540,4510,4480, 4450,4420,4390,4360, 4330,4300,
         4270,4240,4210,4180, 4150,4120,4090,4060,4030,4000,3970,3940,3910,3880,3850,3820,3790,3760,3730,3700,
         3670,3640,3610,3580,3550,3520,3490,3460,3430,3400,3370, 3340, 3310,3280, 3250, 3220, 3190,
         3160, 3130, 3100, 3070, 3040, 3010, 2980, 2950,2920, 2890,2860,2830,2800, 2770,2740,2710, 2680,
         2650, 2620,2590,2560, 2530, 2500,2470,2440,2410, 2380,2350, 2320,2290, 2260, 2230, 2200, 2170, 2140,2110,2080,
         2050,2020,1990,1960,1930,1900,1870,1840,1810,1780,1750,1720,1690,1660,1630,1600,1570,1540,1510,
         1480,1450,1420,1390,1360,1330,1300,1270,1240,1210,1180,1150,1120,1090,1060,1030,1000,970,
         940,910,880,850,820,790,760,730,700,670,640,610,580,550,520,490,460,430,400,370,340,310,280,250, 220)
sums<-cbind(calBP, MKDE, out50)
write.table(sums, file = "data/Sumbin/FremontSumbin.csv", sep = ",", col.names=NA)

##200 year bins
MKDE<-rowMeans(out200)
##Add in the 30 year bin dates
calBP<-c(8000, 7800, 7600, 7400, 7200, 7000, 6800, 6600, 6400, 6200, 6000, 5800, 5600, 5400,
         5200, 5000, 4800, 4600, 4400, 4200, 4000, 3800, 3600, 3400, 3200, 3000, 2800, 2600,
         2400, 2200, 2000, 1800, 1600, 1400, 1200, 1000, 800, 600, 400
)
sums<-cbind(calBP, MKDE, out200)
write.table(sums, file = "data/Sumbin/FremontSumbin200.csv", sep = ",", col.names=NA)


###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

###calculate growth rates for 200 year summed bins

d <-read_csv("data/Sumbin/FremontSumbin200.csv") %>%
  dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))

##Add ID variable for culture history periods
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200, 700, 1850, 4800, 8200),
                         labels=c('Post-Fremont','Fremont','Late Archaic','Archaic'))

write.table(pcgrowth, file = "data/Percapita200/FremontPerCap200.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps.

d <-read.csv("data/Sumbin/FremontSumbin.csv") 
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200, 700, 1850, 4800, 8200),
                         labels=c('Post-Fremont','Fremont','Late Archaic','Archaic'))

write.table(pcgrowth, file = "data/Percapita/FremontPerCap.csv", sep = ",", col.names=NA)


###Plot mean KDE against the per capita growth rate in the North
fre30pc<- read.csv("data/Percapita/FremontPerCap.csv")
#fre30pc<- read.csv("data/Percapita200/FremontPerCap200.csv")


fre30pc2<-subset(fre30pc, calBP<2701 & calBP>400)

StKDE<-(fre30pc2$MKDE-min(fre30pc2$MKDE))/(max((fre30pc2$MKDE)-min(fre30pc2$MKDE)))
fre30pc3<-cbind(StKDE,fre30pc2)

##Write food production file
#write.table(fre30pc3, file = "data/Percapita2/FremontPerCap.csv", sep = ",", col.names=NA)

pcfre <- ggplot(fre30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  scale_y_continuous(limits=c(-.7,0.7))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE density", y="Per capita growth rate", title = "B. Fremont KDE Per Capita Growth vs. Density")+
  geom_hline(yintercept=0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcfre

frecpt <- ggplot(fre30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Time (years cal BP)", y="Standardized KDE density", title = "A. Fremont Density vs. Time")+
  geom_hline(yintercept = 0.20)
#geom_vline(xintercept = 1400)+
#geom_vline(xintercept = 1100)+
#geom_vline(xintercept = 800)+
#annotate("text", x =2000, y = 1, label = "Period 1", size = 6)+
#annotate("text", x =1300, y = 1, label = "2", size = 6)+
#annotate("text", x =950, y = 1, label = "3", size = 6)+
#annotate("text", x =500, y = 1, label = "4", size = 6)
frecpt

Figfremont<-plot_grid(pcfre, frecpt, ncol=2, align="hv", axis = "rl")
Figfremont

#pdf("data/figs/Exfremont.pdf", width=20.55, height=14)
#Figfremont
#dev.off()


#=Chihuahua Desert (Jornada US) Case ID 4===============================================================
#===========================================================
#Load radiocarbon
SPD<-read.csv(file="data/RawP3Kc14.csv", header=T)
#subset the data to the region of study
boxsd<- subset(SPD, Latitude>27 & Latitude<34 & Longitude>-108 & Longitude< -103)
#write.table(boxsd, file = "data/Jornadadates.csv", sep = ",", col.names=NA)
#boxsd<-read.csv(file="data/Jornadadates.csv", header=T)
###MAP
counties<-map_data("state")

ArchGlobeMap<-ggplot() +
  geom_polygon(data = counties, mapping = aes(x = long, y = lat, group = group),
               fill = "grey", color = "white") +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45)+
  #geom_polygon(data = canada, aes(x=long, y = lat, group = group),
  #    fill = "white", color="black") +
  theme_bw()+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Longitude", y="Latitude", title = "US ArchaeoGlobe Regions and Radiocarbon")+
  geom_point(data=boxsd, aes(Longitude, Latitude, color=factor(Province)),
             inherit.aes = FALSE, alpha = 0.5, size = 2)
ArchGlobeMap

boxsd <- boxsd[!is.na(boxsd$SiteID), ]


CalMz <- calibrate(x = boxsd$Age,  errors = boxsd$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = boxsd$SiteID, ages = boxsd$Age, h = 100)

####Run SPD
spd.mz <- spd(CalMz, bins=boxbins, runm=200, timeRange=c(8200,200))
plot(spd.mz, runm=200, xlim=c(8200,200), type="simple")

##Check the effect of h function on SPD if you desire
#binsense(x=CalMz,y=SPD$SiteName,h=seq(0,500,100),timeRange=c(4000,100))

##KDE
####make KDEs
US.randates = sampleDates(CalMz, bins=boxbins, nsim=500,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(8200,200),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')
D.ckde$timeRange

##Write matrix of KDEs as a data frame
Check<-as.data.frame(D.ckde$res.matrix)
#I then convert NAs at the end due to KDE bandwidth to zeros to easily cbind the data frames below
Check2 <- replace(Check, is.na(Check), 0)
#Calculate the mean KDE of each time step of the 200 simulations
MKDE<-rowMeans(Check2)
#calculate the 5th and 95th percentile of each time step of the 200 KDE simulations
lo<-apply( Check2, #select columns
           1, # row-wise calcs
           quantile, probs=0.05) # give `quantile`
hi <-apply( Check2, #select columns
            1, # row-wise calcs
            quantile, probs=0.95) # give `quantile`

calBP<-spd.mz$grid$calBP
PrDens<-spd.mz$grid$PrDens

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<- subset(dd, MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/JornadaKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/JornadaKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')
out200<-rollapply(dd2c,200,(sum),by=200,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
calBP<-c(8170, 8140, 8110, 8080, 8050, 8020, 7990, 7960, 7930, 7900,7870,7840,7810,7780,7750,7720,
         7690,7660,7630,7600,7570,7540,7510,7480,7450,7420,7390,7360,7330,7300,7270,7240,7210,7180, 7150, 7120, 7090,
         7060, 7030, 7000, 6970, 6940, 6910, 6880, 6850, 6820, 6790, 6760, 6730, 6700, 6670, 6640,
         6610, 6580, 6550, 6520, 6490, 6460, 6430, 6400, 6370, 6340, 6310, 6280, 6250, 6220, 6190, 6160, 6130,6100,
         6070,6040,6010,5980,5950,5920,5890,5860,5830,5800,5770,5740,5710,5680,5650,5620,5590,5560,5530,
         5500,5470, 5440,5410, 5380,5350,5320,5290,5260, 5230,5200,5170,5140,5110,5080,5050,5020,4990, 4960,
         4930,4900,4870,4840,4810,4780,4750,4720,4690,4660,4630,4600,4570,4540,4510,4480, 4450,4420,4390,4360, 4330,4300,
         4270,4240,4210,4180, 4150,4120,4090,4060,4030,4000,3970,3940,3910,3880,3850,3820,3790,3760,3730,3700,
         3670,3640,3610,3580,3550,3520,3490,3460,3430,3400,3370, 3340, 3310,3280, 3250, 3220, 3190,
         3160, 3130, 3100, 3070, 3040, 3010, 2980, 2950,2920, 2890,2860,2830,2800, 2770,2740,2710, 2680,
         2650, 2620,2590,2560, 2530, 2500,2470,2440,2410, 2380,2350, 2320,2290, 2260, 2230, 2200, 2170, 2140,2110,2080,
         2050,2020,1990,1960,1930,1900,1870,1840,1810,1780,1750,1720,1690,1660,1630,1600,1570,1540,1510,
         1480,1450,1420,1390,1360,1330,1300,1270,1240,1210,1180,1150,1120,1090,1060,1030,1000,970,
         940,910,880,850,820,790,760,730,700,670,640,610,580,550,520,490,460,430,400,370,340,310,280,250, 220)
sums<-cbind(calBP, MKDE, out50)

write.table(sums, file = "data/Sumbin/JornadaSumbin.csv", sep = ",", col.names=NA)

##200 year bins
MKDE<-rowMeans(out200)
##Add in the 30 year bin dates
calBP<-c(8000, 7800, 7600, 7400, 7200, 7000, 6800, 6600, 6400, 6200, 6000, 5800, 5600, 5400,
         5200, 5000, 4800, 4600, 4400, 4200, 4000, 3800, 3600, 3400, 3200, 3000, 2800, 2600,
         2400, 2200, 2000, 1800, 1600, 1400, 1200, 1000, 800, 600, 400
)
sums<-cbind(calBP, MKDE, out200)
write.table(sums, file = "data/Sumbin/JornadaSumbin200.csv", sep = ",", col.names=NA)


###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

###calculate growth rates for 200 year summed bins

d <-read_csv("data/Sumbin/JornadaSumbin200.csv") %>%
  dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))

##Code cultural historical periods
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200, 650, 1000, 1800, 3800, 6000, 8200),
                         labels=c('Post-Pueblo','Pueblo','Formative','Early Agricultural Period','Middle Archaic','Early Archaic'))

write.table(pcgrowth, file = "data/Percapita200/JornadaPerCap200.csv", sep = ",", col.names=NA)


###Calculate per capita growth rate of 30 year time steps.

d <-read.csv("data/Sumbin/JornadaSumbin.csv") 
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))

##Code cultural historical periods
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200, 650, 1000, 1800, 3800, 6000, 8200),
                         labels=c('Post-Pueblo','Pueblo','Formative','Early Agricultural Period','Middle Archaic','Early Archaic'))

write.table(pcgrowth, file = "data/Percapita/JornadaPerCap.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North
jor30pc<- read.csv("data/Percapita/JornadaPerCap.csv")
#jor30pc<- read.csv("data/Percapita200/JornadaPerCap200.csv")

jor30pc2<-subset(jor30pc, calBP<4000 & calBP>200)

StKDE<-(jor30pc2$MKDE-min(jor30pc2$MKDE))/(max((jor30pc2$MKDE)-min(jor30pc2$MKDE)))
jor30pc3<-cbind(StKDE, jor30pc2)

##Write food production file
#write.table(jor30pc3, file = "data/Percapita2/JornadaPerCap.csv", sep = ",", col.names=NA)

pcjor <- ggplot(jor30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  scale_y_continuous(limits=c(-.3,0.7))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE density", y="KDE per capita growth", title = "A. Jornada KDE Per Capita Growth vs. Density")+
  geom_vline(xintercept = 0.20)+
  geom_hline(yintercept = 0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcjor

jorcpt <- ggplot(jor30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "B. Jornada Density vs. Time")+
  geom_hline(yintercept = 0.20)
#geom_vline(xintercept = 2250)+
#geom_vline(xintercept = 830)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
jorcpt

#paired plot===============================
Figjornada<-plot_grid(pcjor, jorcpt, ncol=2, align="hv", axis = "rl")
Figjornada

#pdf("data/Figs/Exjornada.pdf", width=20.55, height=14)
#Figjornada
#dev.off()

###===Southern Colorado Plateau (Upland US Southwest) Case ID 6===========================

SPD<-read.csv(file="data/RawP3Kc14.csv", header=T)
boxsd<- subset(SPD, Latitude>34 & Latitude<38.01 & Longitude>-113 & Longitude< -105)
#write.table(boxsd, file = "data/ColoradoPlatdates.csv", sep = ",", col.names=NA)
#boxsd<-read.csv(file="data/ColoradoPlatdates.csv", header=T)
###MAP
counties<-map_data("state")

ArchGlobeMap<-ggplot() +
  geom_polygon(data = counties, mapping = aes(x = long, y = lat, group = group),
               fill = "grey", color = "white") +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45)+
  #geom_polygon(data = canada, aes(x=long, y = lat, group = group),
  #    fill = "white", color="black") +
  theme_bw()+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Longitude", y="Latitude", title = "US ArchaeoGlobe Regions and Radiocarbon")+
  geom_point(data=boxsd, aes(Longitude, Latitude, color=factor(Province)),
             inherit.aes = FALSE, alpha = 0.5, size = 2)
ArchGlobeMap

boxsd <- boxsd[!is.na(boxsd$SiteID), ]

CalMz <- calibrate(x = boxsd$Age,  errors = boxsd$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = boxsd$SiteID, ages = boxsd$Age, h = 100)

####Run SPD
spd.mz <- spd(CalMz, bins=boxbins, runm=200, timeRange=c(8200,200))
plot(spd.mz, runm=200, xlim=c(8200,200), type="simple")

calBP<-spd.mz$grid$calBP
PrDens<-spd.mz$grid$PrDens

##Check the effect of h function on SPD if you desire
#binsense(x=CalMz,y=SPD$SiteName,h=seq(0,500,100),timeRange=c(4000,100))

##KDE
####make KDEs
US.randates = sampleDates(CalMz, bins=boxbins, nsim=500,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(8200,200),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')


##Write matrix of KDEs as a data frame
Check<-as.data.frame(D.ckde$res.matrix)
#I then convert NAs at the end due to KDE bandwidth to zeros to easily cbind the data frames below
Check2 <- replace(Check, is.na(Check), 0)
#Calculate the mean KDE of each time step of the 200 simulations
MKDE<-rowMeans(Check2)
#calculate the 5th and 95th percentile of each time step of the 200 KDE simulations
lo<-apply( Check2, #select columns
           1, # row-wise calcs
           quantile, probs=0.05) # give `quantile`
hi <-apply( Check2, #select columns
            1, # row-wise calcs
            quantile, probs=0.95) # give `quantile`

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<- subset(dd, MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/ColoradoPlatKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/ColoradoPlatKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')
out200<-rollapply(dd2c,200,(sum),by=200,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
calBP<-c(8170, 8140, 8110, 8080, 8050, 8020, 7990, 7960, 7930, 7900,7870,7840,7810,7780,7750,7720,
         7690,7660,7630,7600,7570,7540,7510,7480,7450,7420,7390,7360,7330,7300,7270,7240,7210,7180, 7150, 7120, 7090,
         7060, 7030, 7000, 6970, 6940, 6910, 6880, 6850, 6820, 6790, 6760, 6730, 6700, 6670, 6640,
         6610, 6580, 6550, 6520, 6490, 6460, 6430, 6400, 6370, 6340, 6310, 6280, 6250, 6220, 6190, 6160, 6130,6100,
         6070,6040,6010,5980,5950,5920,5890,5860,5830,5800,5770,5740,5710,5680,5650,5620,5590,5560,5530,
         5500,5470, 5440,5410, 5380,5350,5320,5290,5260, 5230,5200,5170,5140,5110,5080,5050,5020,4990, 4960,
         4930,4900,4870,4840,4810,4780,4750,4720,4690,4660,4630,4600,4570,4540,4510,4480, 4450,4420,4390,4360, 4330,4300,
         4270,4240,4210,4180, 4150,4120,4090,4060,4030,4000,3970,3940,3910,3880,3850,3820,3790,3760,3730,3700,
         3670,3640,3610,3580,3550,3520,3490,3460,3430,3400,3370, 3340, 3310,3280, 3250, 3220, 3190,
         3160, 3130, 3100, 3070, 3040, 3010, 2980, 2950,2920, 2890,2860,2830,2800, 2770,2740,2710, 2680,
         2650, 2620,2590,2560, 2530, 2500,2470,2440,2410, 2380,2350, 2320,2290, 2260, 2230, 2200, 2170, 2140,2110,2080,
         2050,2020,1990,1960,1930,1900,1870,1840,1810,1780,1750,1720,1690,1660,1630,1600,1570,1540,1510,
         1480,1450,1420,1390,1360,1330,1300,1270,1240,1210,1180,1150,1120,1090,1060,1030,1000,970,
         940,910,880,850,820,790,760,730,700,670,640,610,580,550,520,490,460,430,400,370,340,310,280,250, 220)
sums<-cbind(calBP, MKDE, out50)
write.table(sums, file = "data/Sumbin/ColoradoPlatSumbin.csv", sep = ",", col.names=NA)

##200 year bins
MKDE<-rowMeans(out200)
##Add in the 30 year bin dates
calBP<-c(8000, 7800, 7600, 7400, 7200, 7000, 6800, 6600, 6400, 6200, 6000, 5800, 5600, 5400,
         5200, 5000, 4800, 4600, 4400, 4200, 4000, 3800, 3600, 3400, 3200, 3000, 2800, 2600,
         2400, 2200, 2000, 1800, 1600, 1400, 1200, 1000, 800, 600, 400
)
sums<-cbind(calBP, MKDE, out200)
write.table(sums, file = "data/Sumbin/ColoradoPlatSumbin200.csv", sep = ",", col.names=NA)
###calculate growth rates for 200 year summed bins

d <-read_csv("data/Sumbin/ColoradoPlatSumbin200.csv") %>%
  dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))

##Code cultural historical periods
##Code cultural historical periods
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200, 320, 670, 820, 1050, 1210, 1455, 3430, 4000, 8200),
                         labels=c('PV','PIV','PIII','PII','PI','BMIII','BMII','Early Agricultural Period','Archaic'))


write.table(pcgrowth, file = "data/Percapita200/ColoradoPlatPerCap200.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps.

d <-read.csv("data/Sumbin/ColoradoPlatSumbin.csv") 
d2<-arrange(d, calBP)
### We calculate per capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))

##Code cultural historical periods
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200, 320, 670, 820, 1050, 1210, 1455, 3430, 4000, 8200),
                         labels=c('PV','PIV','PIII','PII','PI','BMIII','BMII','Early Agricultural Period','Archaic'))

write.table(pcgrowth, file = "data/Percapita/ColoradoPlatPerCap.csv", sep = ",", col.names=NA)


###Plot mean KDE against the per capita growth rate in the North
col30pc<- read.csv("data/Percapita/ColoradoPlatPerCap.csv")
#col30pc<- read.csv("data/Percapita200/ColoradoPlatPerCap200.csv")

col30pc2<-subset(col30pc, calBP<4001 & calBP>200)

StKDE<-(col30pc2$MKDE-min(col30pc2$MKDE))/(max((col30pc2$MKDE)-min(col30pc2$MKDE)))
col30pc3<-cbind(StKDE, col30pc2)

##Write food production file
#write.table(col30pc3, file = "data/Percapita2/ColoradoPlatPerCap.csv", sep = ",", col.names=NA)

pccol <- ggplot(col30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  scale_y_continuous(limits=c(-.3,.7))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Standardized KDE density", y="KDE per capita growth", title = "A. Upland SW US KDE Per Capita Growth vs. Density")+
  geom_hline(yintercept = 0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pccol

colcpt <- ggplot(col30pc3,aes(x=(calBP), y=(StKDE))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=StKDE, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Standardized KDE density", title = "B. Upland SW US Density vs. Time")+
  geom_hline(yintercept = 0.20)
#geom_vline(xintercept = 3320)+
#geom_vline(xintercept = 2550)+
#geom_vline(xintercept = 2250)+
#geom_vline(xintercept = 830)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
colcpt

#paired plot===============================
FiguplandSW<-plot_grid(pccol, colcpt, ncol=2, align="hv", axis = "rl")
FiguplandSW

#pdf("data/Figs/ExuplandSW.pdf", width=20.55, height=14)
#FiguplandSW
#dev.off()

##Load radiocarbon data for Archaeoglobe regional analysis======================================
###Continental US KDE Construction

SPD<-read.csv(file="data/RawNA1.csv", header=T)

SPD <- SPD[!is.na(SPD$SiteID), ]

cptcal <- calibrate(x = SPD$Age,  errors = SPD$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = SPD$SiteID, ages = SPD$Age, h = 100)

####Run analysis for component 3 logistic 3500 to 150
spd.CTx <- spd(cptcal, bins=NA, runm=200, timeRange=c(8200,200))
plot(spd.CTx, runm=200, xlim=c(8200,200), type="simple")

PrDens<-spd.CTx$grid$PrDens
calBP<-spd.CTx$grid$calBP

##Check the effect of h function on SPD if you desire
#binsense(x=CalMz,y=SPD$SiteName,h=seq(0,500,100),timeRange=c(4000,100))

##KDE
####make KDEs
US.randates = sampleDates(cptcal, bins=boxbins, nsim=500,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(8200,200),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')

##Write matrix of KDEs as a data frame
Check<-as.data.frame(D.ckde$res.matrix)
#I then convert NAs at the end due to KDE bandwidth to zeros to easily cbind the data frames below
Check2 <- replace(Check, is.na(Check), 0)
#Calculate the mean KDE of each time step of the 200 simulations
MKDE<-rowMeans(Check2)
#calculate the 10th and 90th percentile of each time step of the 200 KDE simulations
lo<-apply( Check2, #select columns
           1, # row-wise calcs
           quantile, probs=0.1) # give `quantile`
hi <-apply( Check2, #select columns
            1, # row-wise calcs
            quantile, probs=0.9) # give `quantile`

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<-dd %>%  filter(MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/NAKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/NAKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')
out200<- rollapply(dd2c,200,(sum),by=200,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)

#calculate the 5th and 95th percentile of each time step of the 200 KDE simulations
lo<-apply(out50, #select columns
          1, # row-wise calcs
          quantile, probs=0.10) # give `quantile`
hi <-apply(out50, #select columns
           1, # row-wise calcs
           quantile, probs=0.90) # give `quantile`
calBP<-c(8170, 8140, 8110, 8080, 8050, 8020, 7990, 7960, 7930, 7900,7870,7840,7810,7780,7750,7720,
         7690,7660,7630,7600,7570,7540,7510,7480,7450,7420,7390,7360,7330,7300,7270,7240,7210,7180, 7150, 7120, 7090,
         7060, 7030, 7000, 6970, 6940, 6910, 6880, 6850, 6820, 6790, 6760, 6730, 6700, 6670, 6640,
         6610, 6580, 6550, 6520, 6490, 6460, 6430, 6400, 6370, 6340, 6310, 6280, 6250, 6220, 6190, 6160, 6130,6100,
         6070,6040,6010,5980,5950,5920,5890,5860,5830,5800,5770,5740,5710,5680,5650,5620,5590,5560,5530,
         5500,5470, 5440,5410, 5380,5350,5320,5290,5260, 5230,5200,5170,5140,5110,5080,5050,5020,4990, 4960,
         4930,4900,4870,4840,4810,4780,4750,4720,4690,4660,4630,4600,4570,4540,4510,4480, 4450,4420,4390,4360, 4330,4300,
         4270,4240,4210,4180, 4150,4120,4090,4060,4030,4000,3970,3940,3910,3880,3850,3820,3790,3760,3730,3700,
         3670,3640,3610,3580,3550,3520,3490,3460,3430,3400,3370, 3340, 3310,3280, 3250, 3220, 3190,
         3160, 3130, 3100, 3070, 3040, 3010, 2980, 2950,2920, 2890,2860,2830,2800, 2770,2740,2710, 2680,
         2650, 2620,2590,2560, 2530, 2500,2470,2440,2410, 2380,2350, 2320,2290, 2260, 2230, 2200, 2170, 2140,2110,2080,
         2050,2020,1990,1960,1930,1900,1870,1840,1810,1780,1750,1720,1690,1660,1630,1600,1570,1540,1510,
         1480,1450,1420,1390,1360,1330,1300,1270,1240,1210,1180,1150,1120,1090,1060,1030,1000,970,
         940,910,880,850,820,790,760,730,700,670,640,610,580,550,520,490,460,430,400,370,340,310,280,250, 220)
sums<-cbind(calBP, MKDE, hi, lo, out50)

write.table(sums, file = "data/Sumbin/NASumbin.csv", sep = ",", col.names=NA)

##200 year bins
MKDE<-rowMeans(out200)
##Add in the 30 year bin dates
calBP<-c(8000, 7800, 7600, 7400, 7200, 7000, 6800, 6600, 6400, 6200, 6000, 5800, 5600, 5400,
         5200, 5000, 4800, 4600, 4400, 4200, 4000, 3800, 3600, 3400, 3200, 3000, 2800, 2600,
         2400, 2200, 2000, 1800, 1600, 1400, 1200, 1000, 800, 600, 400
)
sums<-cbind(calBP, MKDE, out200)
write.table(sums, file = "data/Sumbin/NASumbin200.csv", sep = ",", col.names=NA)


###calculate growth rates for 200 year summed bins

d <-read_csv("data/Sumbin/NASumbin200.csv") %>%
  dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))

pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200, 500, 790, 2000, 4200, 8200),
                         labels=c('Historic','Transition','Agriculture','Early Agriculture','Archaic'))


write.table(pcgrowth, file = "data/Percapita200/NAPerCap200.csv", sep = ",", col.names=NA)


###Calculate per capita growth rate of 30 year time steps.

d <-read_csv("data/Sumbin/NASumbin.csv") %>%
  dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200, 500, 790, 2000, 4200, 8200),
                         labels=c('Historic','Transition','Agriculture','Early Agriculture','Archaic'))


write.table(pcgrowth, file = "data/Percapita/NAPerCap.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North
ct30pc<- read.csv("data/Percapita/NAPerCap.csv")

ct30pc2<-subset(ct30pc, calBP<8000 & calBP>200)

#Standardize the mean KDE by the maximum mean KDE during the Neolithic 
StKDE<-(ct30pc2$MKDE-min(ct30pc2$MKDE))/(max((ct30pc2$MKDE)-min(ct30pc2$MKDE)))
##Add the standardized KDE to the Neolithic dataframe
ct30pc3<-cbind(StKDE,ct30pc2)

##Write food production file
#write.table(ct30pc3, file = "data/Recession/NorthAmericaPercapRecession.csv", sep = ",", col.names=NA)

pcctex <- ggplot(ct30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  #scale_y_continuous(limits=c(-.1,0.25))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Mean KDE (density)", y="KDE per capita growth", title = "A. Continuous Per Capita Growth vs. Density")+
  geom_hline(yintercept=0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcctex

ctexcpt <- ggplot(ct30pc3,aes(x=(calBP), y=(MKDE*100))) +
  geom_ribbon(aes(ymin = lo*100, ymax = hi*100), fill = "grey70") +
  geom_point(aes(y=MKDE*100, color=factor(PeriodID)), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Mean KDE density", title = "B. Continuous US Mean KDE vs. Time")
#geom_hline(yintercept=0.25)
ctexcpt

#paired plot===============================
Figctex<-plot_grid(pcctex, ctexcpt, ncol=2, align="hv", axis = "rl")
Figctex


#SPlit data for SW region============================================
SPD.1<- subset(SPD, Archaeo_ID=="15")

cptcal <- calibrate(x = SPD.1$Age,  errors = SPD.1$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = SPD.1$SiteID, ages = SPD.1$Age, h = 100)

####Run analysis for component 3 logistic 3500 to 150
spd.CTx <- spd(cptcal, bins=NA, runm=200, timeRange=c(8200,200))
plot(spd.CTx, runm=200, xlim=c(8200,200), type="simple")

PrDens<-spd.CTx$grid$PrDens
calBP<-spd.CTx$grid$calBP

##Check the effect of h function on SPD if you desire
#binsense(x=CalMz,y=SPD$SiteName,h=seq(0,500,100),timeRange=c(4000,100))

##KDE
####make KDEs
US.randates = sampleDates(cptcal, bins=boxbins, nsim=500,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(8200,200),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')

##Write matrix of KDEs as a data frame
Check<-as.data.frame(D.ckde$res.matrix)
#I then convert NAs at the end due to KDE bandwidth to zeros to easily cbind the data frames below
Check2 <- replace(Check, is.na(Check), 0)
#Calculate the mean KDE of each time step of the 200 simulations
MKDE<-rowMeans(Check2)
#calculate the 10th and 90th percentile of each time step of the 200 KDE simulations
lo<-apply( Check2, #select columns
           1, # row-wise calcs
           quantile, probs=0.1) # give `quantile`
hi <-apply( Check2, #select columns
            1, # row-wise calcs
            quantile, probs=0.9) # give `quantile`

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<-dd %>%  filter(MKDE >0)
##Write the table
write.table(dd2, file = "data/KDEs/SW15KDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("data/KDEs/SW15KDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')
out200<- rollapply(dd2c,200,(sum),by=200,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)

#calculate the 5th and 95th percentile of each time step of the 200 KDE simulations
lo<-apply(out50, #select columns
          1, # row-wise calcs
          quantile, probs=0.10) # give `quantile`
hi <-apply(out50, #select columns
           1, # row-wise calcs
           quantile, probs=0.90) # give `quantile`
calBP<-c(8170, 8140, 8110, 8080, 8050, 8020, 7990, 7960, 7930, 7900,7870,7840,7810,7780,7750,7720,
         7690,7660,7630,7600,7570,7540,7510,7480,7450,7420,7390,7360,7330,7300,7270,7240,7210,7180, 7150, 7120, 7090,
         7060, 7030, 7000, 6970, 6940, 6910, 6880, 6850, 6820, 6790, 6760, 6730, 6700, 6670, 6640,
         6610, 6580, 6550, 6520, 6490, 6460, 6430, 6400, 6370, 6340, 6310, 6280, 6250, 6220, 6190, 6160, 6130,6100,
         6070,6040,6010,5980,5950,5920,5890,5860,5830,5800,5770,5740,5710,5680,5650,5620,5590,5560,5530,
         5500,5470, 5440,5410, 5380,5350,5320,5290,5260, 5230,5200,5170,5140,5110,5080,5050,5020,4990, 4960,
         4930,4900,4870,4840,4810,4780,4750,4720,4690,4660,4630,4600,4570,4540,4510,4480, 4450,4420,4390,4360, 4330,4300,
         4270,4240,4210,4180, 4150,4120,4090,4060,4030,4000,3970,3940,3910,3880,3850,3820,3790,3760,3730,3700,
         3670,3640,3610,3580,3550,3520,3490,3460,3430,3400,3370, 3340, 3310,3280, 3250, 3220, 3190,
         3160, 3130, 3100, 3070, 3040, 3010, 2980, 2950,2920, 2890,2860,2830,2800, 2770,2740,2710, 2680,
         2650, 2620,2590,2560, 2530, 2500,2470,2440,2410, 2380,2350, 2320,2290, 2260, 2230, 2200, 2170, 2140,2110,2080,
         2050,2020,1990,1960,1930,1900,1870,1840,1810,1780,1750,1720,1690,1660,1630,1600,1570,1540,1510,
         1480,1450,1420,1390,1360,1330,1300,1270,1240,1210,1180,1150,1120,1090,1060,1030,1000,970,
         940,910,880,850,820,790,760,730,700,670,640,610,580,550,520,490,460,430,400,370,340,310,280,250, 220)
sums<-cbind(calBP, MKDE, hi, lo, out50)

write.table(sums, file = "data/Sumbin/SW15Sumbin.csv", sep = ",", col.names=NA)

##200 year bins
MKDE<-rowMeans(out200)
##Add in the 30 year bin dates
calBP<-c(8000, 7800, 7600, 7400, 7200, 7000, 6800, 6600, 6400, 6200, 6000, 5800, 5600, 5400,
         5200, 5000, 4800, 4600, 4400, 4200, 4000, 3800, 3600, 3400, 3200, 3000, 2800, 2600,
         2400, 2200, 2000, 1800, 1600, 1400, 1200, 1000, 800, 600, 400
)
sums<-cbind(calBP, MKDE, out200)
write.table(sums, file = "data/Sumbin/SW15Sumbin200.csv", sep = ",", col.names=NA)


###calculate growth rates for 200 year summed bins

d <-read_csv("data/Sumbin/SW15Sumbin200.csv") %>%
  dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))

pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200, 500, 790, 2000, 4200, 8200),
                         labels=c('Historic','Transition','Agriculture','Early Agriculture','Archaic'))


write.table(pcgrowth, file = "data/Percapita200/SW15PerCap200.csv", sep = ",", col.names=NA)


###Calculate per capita growth rate of 30 year time steps.

d <-read_csv("data/Sumbin/SW15Sumbin.csv") %>%
  dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))
pcgrowth$PeriodID <- cut(pcgrowth$calBP,
                         breaks=c(200, 500, 790, 2000, 4200, 8200),
                         labels=c('Historic','Transition','Agriculture','Early Agriculture','Archaic'))


write.table(pcgrowth, file = "data/Percapita/SW15PerCap.csv", sep = ",", col.names=NA)

###Plot mean KDE against the per capita growth rate in the North
ct30pc<- read.csv("data/Percapita/SW15PerCap.csv")
ct30pc2<-subset(ct30pc, calBP<4000 & calBP>200)

#Standardize the mean KDE by the maximum mean KDE during the Neolithic 
StKDE<-(ct30pc2$MKDE-min(ct30pc2$MKDE))/(max((ct30pc2$MKDE)-min(ct30pc2$MKDE)))
##Add the standardized KDE to the Neolithic dataframe
ct30pc3<-cbind(StKDE,ct30pc2)

##Write food production file
#write.table(ct30pc3, file = "data//Southwest_15_PercapRecession.csv", sep = ",", col.names=NA)

pcctex <- ggplot(ct30pc3,aes(x=(StKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=PeriodID), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  #scale_y_continuous(limits=c(-.1,0.25))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Mean KDE (density)", y="KDE per capita growth", title = "A. American Southwest Per Capita Growth vs. Density")+
  geom_hline(yintercept=0)
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcctex

ctexcpt <- ggplot(ct30pc3,aes(x=(calBP), y=(MKDE*100))) +
  geom_ribbon(aes(ymin = lo*100, ymax = hi*100), fill = "grey70") +
  geom_point(aes(y=MKDE*100, color=factor(PeriodID)), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse()+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years cal BP", y="Mean KDE density", title = "B. American Southwest Mean KDE vs. Time")
#geom_hline(yintercept=0.25)
ctexcpt

#paired plot===============================
FigSW<-plot_grid(pcctex, ctexcpt, ncol=2, align="hv", axis = "rl")
FigSW

