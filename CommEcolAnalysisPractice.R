# Comm Ecol Analysis Practice
library(picante)


#### Indicator Species Analysis #####
library(indicspecies)

# Example from : https://cran.r-project.org/web/packages/indicspecies/vignettes/indicspeciesTutorial.pdf
data(wetland)
# Sites are rows, species are columns
# Create groups on which to partition data
groups = c(rep(1, 17), rep(2, 14), rep(3,10))
# Calculate indicator species associations on wetland species/sites based on the 3 groups
# how(nperm) does a bunch of permutations of this function so a p-value can be calculated
indval = multipatt(wetland, groups,control = how(nperm=999))
summary(indval)
# Indicator species analysis is conducted for each species independently -- therefore, need to
# KEEP IN MIND ISSUES WITH MULTIPLE TESTING (if you test a bunch of sp. doing this, 5% expected to come back positive)
# Maureen berg did a Sidak correction? How do I do this?

# Look at the different components (A and B) of the association index
# A tells you the likelihood of belonging to a certain group (clutch) given that this species is present
# B tells you how likely you are to find this species in a certain group (clutch)
summary(indval, indvalcomp=TRUE)

# If you want to display indicator value summary for all species (not just significant ones):
summary(indval, alpha = 1)
# Otherwise the default alpha = 0.05
# This will still hide species that occur in all groups though (since there's no control group to compare them to - no way to perform statistical test)
# to look at this, use:
indval$sign  # and see "NA" for pvalue
