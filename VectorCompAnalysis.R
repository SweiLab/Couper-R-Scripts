# Vector competence project analysis 
# Part 1: bayesian credible intervals for acquisition experiment
# Part 2: Time course for tick attachment procedure

#### PART 1: ACQUISITION EXPERIMENT ####
library(BEST)
library(BayesianFirstAid)

# Will be using the bayes.prop.test function from the Bayesian First Aid package
# Using example from here: http://www.sumsar.net/blog/2014/06/bayesian-first-aid-prop-test/

# For acquisition comparison:
# Ipacs had 22/30 positives for CA4 and 6/30 positives for B31
# Iscaps had 6/30 positives for CA4 and 22/30 positives for B31

# First doing frequentist approach
Ipacn_pos = c(22,6)
Ipacn_samples = c(30,30)
Ipacs = prop.test(Ipacn_pos, Ipacn_samples)
IpacsBayesian = bayes.prop.test(Ipacn_pos, Ipacn_samples)
# 95% credibe interval for CA4 for Ipacs: [0.56 0.86]
# 95$ credible interval for B31 for Ipacs: [0.085, 0.36]

Iscaps_pos = c(6,22)
Iscaps_samples = c(30,30)
Iscaps = prop.test(Iscaps_pos, Iscaps_samples)
IscapsBayesian = bayes.prop.test(Iscaps_pos, Iscaps_samples)
# Same credible intervals calcualted (which makes sense)

# Plot with credible intervals
#density1 <- seq(5,35,length.out=7)
Values = c(73.3, 20, 20, 73.3)
MyPlot = barplot(Values, ylim = c(0,100), 
                 col  = c("chartreuse4", "blue4", "chartreuse4", "blue4"),
                 space = c(0.1, 0.07, 0.5, 0.07), cex.axis = 1.7)
                 #density = c(20,100,20,100))
Mins = c(56, 8.6, 8.6, 56)
Maxs = c(86, 36, 36, 85)
segments(MyPlot, Mins, MyPlot, Maxs, lwd = 1.5)
arrows(MyPlot, Mins, MyPlot, Maxs, lwd = 1.5, code =3, angle = 90, length = 0.05)

#### PART 2: ATTACHMENT COMPARISON ######
atc = read.csv("~/Documents/Vector_Comp_Project/VectorComp_AttachmentTimeCourse.csv", header= TRUE)
rownames(atc) = atc[,1]
atc= atc[,-1]
hours = c(24,48,72,96,120,144)
days = c(1,2,3,4,5,6)

plot(y= atc[,1], x=days, type = "o",col = "green3", xlab = "Days", ylab = "Drop-offs", lwd = 3,
     main = "Attachment Rates", cex.axis = 1.5)
lines(atc[,2], x= days, type = "o", col = "green4", lwd = 3)
lines(atc[,3], x = days, type = "o", col = "lightskyblue3", lwd = 3)
lines(atc[,4], x = days, type = "o", col = "mediumblue", lwd = 3)


atc2 = read.csv("~/Documents/Vector_Comp_Project/VectorComp_Attachment_CD.csv", header= TRUE)
rownames(atc2) = atc2[,1]
atc2= atc2[,-1]
days = c(1,2,3,4,5,6)

plot(y= atc2[,1], x=days, type = "o",col = "green3", xlab = "Days", ylab = "Drop-offs", lwd = 3,
     main = "Attachment Rates", cex.axis = 1.5)
lines(atc2[,2], x= days, type = "o", col = "green4", lwd = 3)
lines(atc2[,3], x = days, type = "o", col = "lightskyblue3", lwd = 3)
lines(atc2[,4], x = days, type = "o", col = "mediumblue", lwd = 3)


# Using prop test to evaluate statistical significance of this difference
# A chi-squared test?
dropoffs = c(102,20)
placed = c(300,300)
prop.test(dropoffs, placed, p = NULL)


# Weights comparison
wts = read.csv("~/Documents/Vector_Comp_Project/18VectorComp_TickWeights.csv", header = TRUE) 
