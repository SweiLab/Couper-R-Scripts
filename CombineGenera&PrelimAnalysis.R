# Genera Combine Script, Barplot of Genera by Time, and Community Ecology Calculations

library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(randomcoloR)
library(vegan)


### Read in File, reshape/reformat it and create boxplots by time

pwc = read.csv("~/Documents/Pepperwood_Project/Data_Files/Pepperwood_RareRemoved_OTUs_decontam.csv", header = TRUE)
pwc = aggregate(pwc[,2:69], by = list(pwc$Genera), FUN= sum)
# Now only 23 genera remain 

for (i in 2:69)
{pwc[,i] = as.numeric(as.character(pwc[,i]))}

pwc2 = pwc  
# Need to convert read counts into proportions
for (j in 2:69)
{for (i in 1:23)
{pwc2[i,j] = (pwc[i,j]/10182)*100}}

pwc = pwc2
pcprops = pwc

for (j in 2:69)
  for (i in 1:23)
    pcprops[i,j] = (pwc[i,j]/100)


# Props with no adults
pcprops_noadults = pcprops[, -c(2,22,46)]
range(pcprops_noadults[19,-1])
      

      
##### Boxplots showing genera over time #####

Clutch = c("Adult1", rep("1", 19), "Adult2", rep("2", 23), "Adult3", rep("3", 23))
Time = c("Adults", 0,  rep("2", 5), rep("4", 6), rep("6", 7),
         "Adults", rep("0", 6), rep("2", 7), rep("4", 7), rep("6", 3), 
         "Adults", rep("0", 5), 2, 2, rep("4", 8), rep("6", 8))


pwct = t(pwc)
colnames(pwct) = pwct[1,]
pwct = pwct[-1,]
pwct = cbind(pwct, Time)
pwct = as.data.frame(pwct)

# Change class to numeric before aggregating
for (i in 1:23)
{pwct[,i] = as.numeric(as.character(pwct[,i]))}

# Want to sum counts for each genera based on time
pwcta = aggregate(pwct[,1:23], by = list(pwct$Time), FUN = mean)

# Now transpose back to original form
pwcgo = t(pwcta)
colnames(pwcgo) = pwcgo[1,]
columnnames = colnames(pwcgo)
pwcgo = pwcgo[-1,]
pwcgo = cbind(rownames(pwcgo), pwcgo)
rownames(pwcgo) = c()
colnames(pwcgo) = c("Genus", "T0", "T2", "T4", "T6", "Adults")

pwcgo = as.data.frame(pwcgo)
for (i in 2:6)
{pwcgo[,i] = as.numeric(as.character(pwcgo[,i]))}

# Now want to convert between wide and long form
pwcm = melt(pwcgo, id = c("Genus"))
colnames(pwcm) = c("Genera", "Time", "Reads")

# Now plot
# Manually picking colors to make them most distinct

mypal2 = colorRampPalette( brewer.pal(12 , "Paired" ) )
colors= mypal2(23)

#pwm$Time = factor(pwm$Time,levels = c("Adults", "T0", "T2", "T4", "T6"))
pwcm =  pwcm[with(pwcm, order(pwcm$Time, -pwcm$Reads)), ]

#pwcm$Genera <- factor(pwcm$Genera, levels = pwcm$Genera)

p = ggplot(pwcm, aes(x= pwcm$Time, y=pwcm$Reads, fill= reorder(pwcm$Genera, pwcm$Reads)))+
  geom_bar(stat="identity", colour='black')+        # contour colour
  scale_fill_manual(values = sample(colors,23) )+ 
  xlab("Environmental Exposure")+                                     
  ylab("Proportion of Reads") +
  theme(legend.text=element_text(size=18))
p + guides(fill=guide_legend(reverse=TRUE), col= guide_legend(nrow=1)) 

# See Pepperwood Manuscript Plots R script for final figure


#### Comm Ecol Calculations: Species Richness ####
# Go back to original dataset for these calculations
SpeciesRich = rep(0,68)
for (i in 2:69)
{SpeciesRich[i-1] = sum(pwc[,i] > 0)}
Clutch = c("Adult1", rep("1", 19), "Adult2", rep("2", 23), "Adult3", rep("3", 23))
Time = c("Adult1", 0,  rep("2", 5), rep("4", 6), rep("6", 7),
         "Adult2", rep("0", 6), rep("2", 7), rep("4", 7), rep("6", 3), 
         "Adult3", rep("0", 5), 2, 2, rep("4", 8), rep("6", 8))
SpeciesRichClutch = cbind(SpeciesRich, Clutch)
SpeciesRichClutch = as.data.frame(SpeciesRichClutch)
SpeciesRichClutch$SpeciesRich = as.numeric(as.character(SpeciesRichClutch$SpeciesRich))
SpRichAvgs = aggregate(SpeciesRichClutch$SpeciesRich, by = list(SpeciesRichClutch$Clutch), FUN = mean) 
colnames(SpRichAvgs) = c("Clutch", "Mean")
SpRichSd = aggregate(SpeciesRichClutch$SpeciesRich, by = list(SpeciesRichClutch$Clutch), FUN = sd)
colnames(SpRichSd) = c("Clutch", "Sd")
SpRichCSum = cbind(SpRichAvgs, SpRichSd$Sd)
colnames(SpRichCSum) = c("Clutch", "Mean", "Sd")
# Remove Adults
srcs = SpRichCSum[-c(4:6),]
# Run ANOVA
# Remove adults from original
sc_noadults = SpeciesRichClutch[-c(1,21,45),]
fit = aov(SpeciesRich ~ Clutch, data = sc_noadults)
summary(fit)
# Tukeys post-hoc test (To tell you which clutches actually differ)
posthoc <- TukeyHSD(x=fit, 'Clutch', conf.level=0.95)
# Result = significant difference in species richness between clutches (P<0.001)

SpeciesRichTime = cbind(SpeciesRich, Time)
SpeciesRichTime = as.data.frame(SpeciesRichTime)
SpeciesRichTime$SpeciesRich = as.numeric(as.character(SpeciesRichTime$SpeciesRich))
SpRichTAvgs = aggregate(SpeciesRichTime$SpeciesRich, by = list(SpeciesRichTime$Time), FUN = mean) 
colnames(SpRichTAvgs) = c("Time", "Mean")
SpRichTSd = aggregate(SpeciesRichTime$SpeciesRich, by = list(SpeciesRichTime$Time), FUN = sd)
colnames(SpRichTSd) = c("Time", "Sd")
SpRichTSum = cbind(SpRichTAvgs, SpRichTSd$Sd)
colnames(SpRichTSum) = c("Time", "Mean", "Sd")
# Remove Adults
srts = SpRichTSum[-c(5:7),]
# Run ANOVA
# Remove adults from original
st_noadults = SpeciesRichTime[-c(1,21,45),]
fit2 = aov(SpeciesRich ~ Time, data = st_noadults)
summary(fit2)
posthoc <- TukeyHSD(x=fit2, 'Time', conf.level=0.95)

##### Community Ecology Calculations #####


# Result = significant difference in species richness by time points (P= 0.0117)

# SpRich Between Adults & Larvae
AgeStatus = c("Adult", rep("Larvae", 19), "Adult", rep("Larvae",23), "Adult", rep("Larvae",23))
sr_age = cbind(SpeciesRichClutch, AgeStatus)
sr_aggX = aggregate(sr_age$SpeciesRich, by = list(sr_age$AgeStatus), FUN = mean)
sr_aggY = aggregate(sr_age$SpeciesRich, by = list(sr_age$AgeStatus), FUN = sd)
srSum = cbind(sr_aggX, sr_aggY$x)
colnames(srSum) = c("LifeStage","Mean", "Sd")
# Run ANOVA
fit3 = aov(SpeciesRich ~ AgeStatus, data = sr_age)


#### Comm Ecol Calculations: Species Diversity ####

Divcp = pcprops

#for (j in 2:69)
#{for (i in 1:27)
#{Divp[i,j] = (pw2[i,j]/100)}}

Divcp3 = Divcp
for (j in 2:69)
for (i in 1:23)
  Divcp3[i,j] = (pcprops[i,j]*log(pcprops[i,j]))

ShannonsD = rep(0,68)
for (i in 2:69)
{ShannonsD[i-1] = -1*sum(Divcp3[,i], na.rm = TRUE)}

# Species Diversity by Clutch #
DivClutch = cbind(ShannonsD, Clutch)
DivClutch = as.data.frame(DivClutch)
DivClutch$ShannonsD = as.numeric(as.character(DivClutch$ShannonsD))
dccagg = aggregate(DivClutch$ShannonsD, by = list(DivClutch$Clutch), FUN = mean)
dccaggsd = aggregate(DivClutch$ShannonsD, by = list(DivClutch$Clutch), FUN = sd)
DCcsum = cbind(dccagg, dccaggsd$x)
colnames(DCcsum) = c("Clutch", "Mean", "Sd")
# Remove Adults from original for ANOVA
dcc2 = DivClutch[-c(1,21,45),]
fit4 = aov(ShannonsD ~ Clutch, data = dcc2)
summary(fit4)
# Result = significant differences by clutch (P<0.001)
posthoc2 <- TukeyHSD(x=fit4, 'Clutch', conf.level=0.95)

# Species Diversity by Time #
DivTime = cbind(ShannonsD, Time)
DivTime = as.data.frame(DivTime)
DivTime$ShannonsD = as.numeric(as.character(DivTime$ShannonsD))
dtcagg = aggregate(DivTime$ShannonsD, by = list(DivTime$Time), FUN = mean)
dtcaggsd = aggregate(DivTime$ShannonsD, by = list(DivTime$Time), FUN = sd)
DTcsum = cbind(dtcagg, dtcaggsd$x)
colnames(DTcsum) = c("Time", "Mean", "Sd")
# Remove Adults from original for ANOVA
dtc2 = DivTime[-c(1,21,45),]
fit5 = aov(ShannonsD ~ Time, data = dtc2)
summary(fit5)
posthoc3 <- TukeyHSD(x=fit5, 'Time', conf.level=0.95)
# Result = significant differences by time (P = 0.0108)
boxplot(ShannonsD ~ Time, data = dtc2, main = "Diversity Over Time")

# Species Diversity by life stage (adult/larvae) #
AgeStatus = c("Adult", rep("Larvae", 19), "Adult", rep("Larvae",23), "Adult", rep("Larvae",23))
sd_age = cbind(DivClutch, AgeStatus)
sd_aggX = aggregate(sd_age$ShannonsD, by = list(sd_age$AgeStatus), FUN = mean)
sd_aggY = aggregate(sd_age$ShannonsD, by = list(sd_age$AgeStatus), FUN = sd)
sdcSum = cbind(sd_aggX, sd_aggY$x)
colnames(sdcSum) = c("LifeStage","Mean", "Sd")
# Run ANOVA
fit6 = aov(ShannonsD ~ AgeStatus, data = sd_age)

#### Comm Ecol Calculations: Species Evenness Calculations ####
# Evenness = Shannons Diversity / ln(Species Richness)
Even = rep(0,68)

for (i in 1:68)
{Even[i] = DivClutch$ShannonsD[i]/log(SpeciesRichClutch$SpeciesRich[i])}

# Evenness by Clutch
Even = as.numeric(as.character(Even))
EvenClutch = cbind(Even, Clutch)
EvenClutch = as.data.frame(EvenClutch)
EvenClutch$Even = as.numeric(as.character(EvenClutch$Even))
ECagg = aggregate(EvenClutch$Even, by = list(EvenClutch$Clutch), FUN = mean)
ECaggsd = aggregate(EvenClutch$Even, by = list(EvenClutch$Clutch), FUN = sd)
ECsum = cbind(ECagg, ECaggsd$x)
colnames(ECsum) = c("Clutch", "Mean", "Sd")
# Remove Adults
EC2 = EvenClutch[-c(1,21,45),]
fit7 = aov(Even ~ Clutch, data = EC2)
summary(fit7)
posthoc4 <- TukeyHSD(x=fit7, 'Clutch', conf.level=0.95)

# Evenness by Time
Even = as.numeric(as.character(Even))
EvenTime = cbind(Even, Time)
EvenTime = as.data.frame(EvenTime)
EvenTime$Even = as.numeric(as.character(EvenTime$Even))
ETagg = aggregate(EvenTime$Even, by = list(EvenTime$Time), FUN = mean)
ETaggsd = aggregate(EvenTime$Even, by = list(EvenTime$Time), FUN = sd)
ETsum = cbind(ETagg, ETaggsd$x)
colnames(ETsum) = c("Time", "Mean", "Sd")
# Remove Adults
ET2 = EvenTime[-c(1,21,45),]
fit8 = aov(Even ~ Time, data = ET2)
summary(fit8)

# Evenness by Life Stage
AgeStatus = c("Adult", rep("Larvae", 19), "Adult", rep("Larvae",23), "Adult", rep("Larvae",23))
e_age = cbind(EvenClutch, AgeStatus)
e_aggX = aggregate(e_age$Even, by = list(e_age$AgeStatus), FUN = mean)
e_aggY = aggregate(e_age$Even, by = list(e_age$AgeStatus), FUN = sd)
eSum = cbind(e_aggX, e_aggY$x)
colnames(eSum) = c("LifeStage","Mean", "Sd")
# Run ANOVA
fit9 = aov(Even ~ AgeStatus, data = e_age)
summary(fit9)


# Rickettsia Abundance By Time
Ricks = read.csv("~/Documents/Pepperwood_RickettsiaqPCR.csv", header= TRUE)
# Remove weird extra rows that were added
fit29 = aov(Avg ~ TimePoint, data = Ricks)
summary(fit29)
# No significant difference in Rickettsia over time (despite increase in relative abundance)
boxplot(Avg ~ TimePoint, data = Ricks, main = "Rickettsia Over Time")

# Overall 16s Abundance By Time
Total = read.csv("~/Documents/Pepperwood_16sLoads.csv", header = TRUE)
fit30 = aov(Avg ~ TimePoint, data = Total)
summary(fit30)
boxplot(Avg ~ TimePoint, data = Total, main = "Total Loads Over Time")


# 
