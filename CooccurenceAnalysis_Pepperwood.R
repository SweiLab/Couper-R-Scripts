# Correlation between OTUs using Spearman's rank
# For network inference, we calculated all possible Spearman's rank correlations between OTUs
# Following Barberan 2011 paper

# Example - using ppct data frame created in "CommEcolAnalysisPractice" script
# ppct

library(lessR)
library(igraph)

X1 = ppct[,1] # Acinetobacter
X2 = ppct[,2] # Actinomyces
test = Correlation(X1, X2, method = "spearman", exact = FALSE) # rho = -0.055, p = 0.688
test$pvalue # how to output just the pvalue
test$r # how to output just the rho value
# Barberan paper counted it as significant is >0.6 and p<0.01 

# Now need to scale this up to do all calculations
# First step = remove the Other column
ppctgo =subset(ppct, select=-Other)
# Now have 22 OTUs total - how many pairs is this? (22 choose 2) = 231
# Store rho values in one vector, p-values in another
Rhos = rep(0,231)
Ps = rep(0,231)
# Create indexes from which to pull (to get all combos of OTUs)
indices = combn(22,2)


for (i in 1:231)
{a = indices[1,i]
b = indices[2,i]
a2 = ppctgo[,a]
b2 = ppctgo[,b]
test = Correlation(a2,b2, method = "spearman", exact = FALSE)
Rhos[i] = test$r
Ps[i] = test$pvalue}

# Adjust p-values for multiple testing
Ps2 = p.adjust(Ps, method = "fdr")

# Rhos now contains all the spearman rho values
# Ps now contains all the pvalues

# Look for those which have p< 0.01

# Example of creating a network graph from correlation coefficients 

# library
library(igraph)
# data
head(mtcars)
# Make a correlation matrix:
mat=cor(t(mtcars[,c(1,3:6)]))
# Keep only high correlations
mat[mat<0.995]=0
# Make an Igraph object from this matrix:
network=graph_from_adjacency_matrix( mat, weighted=T, mode="undirected", diag=F)
plot(network)
library(RColorBrewer)
coul = brewer.pal(nlevels(as.factor(mtcars$cyl)), "Set2")
# Map the color to cylinders
my_color=coul[as.numeric(as.factor(mtcars$cyl))]
# plot
par(bg="grey13", mar=c(0,0,0,0))
set.seed(4)
plot(network, 
     vertex.size=12,
     vertex.color=my_color, 
     vertex.label.cex=0.7,
     vertex.label.color="white",
     vertex.frame.color="transparent"
)
# title and legend
text(0,0,"The network chart of the mtcars dataset",col="white")
text(0.2,-0.1," - by the R graph gallery",col="white")
legend(x=-0.6, y=-0.12, legend=paste( levels(as.factor(mtcars$cyl)), " cylinders", sep=""), col = coul , bty = "n", pch=20 , pt.cex = 2, cex = 1, text.col="white" , horiz = T)


# Back to my data - need to populate a data frame with correlation values corresponding to their original OTUs
# First create empty data frame with OTU names from ppct
Corrs = matrix(nrow = 22, ncol = 22)
Corrs = as.data.frame(Corrs)
colnames(Corrs) = colnames(ppctgo)
rownames(Corrs) = colnames(ppctgo)

for (i in 1:231)
{f = indices[1,i]
g = indices[2,i]
Corrs[f,g] = Rhos[i]
Corrs[g,f] = Rhos[i]
Corrs[f,f]= 1}
Corrs[22,22] = 1
Corrs = as.matrix(Corrs)

# Keep only the Rho values above 0.6. If they're less - coerce to 0
Corrs[Corrs<0.6]=0
colnames(Corrs)[1] = "Acinetobacter"

# Now want to do same thing with p-values
Pvals = matrix(nrow = 22, ncol = 22)
Pvals = as.data.frame(Pvals)
colnames(ppctgo)[1] = "Acinetobacter"
colnames(Pvals) = colnames(ppctgo)

rownames(Pvals) = colnames(ppctgo)

for (i in 1:231)
{f = indices[1,i]
g = indices[2,i]
Pvals[f,g] = Ps2[i]
Pvals[g,f] = Ps2[i]
Pvals[f,f]= 0}
Pvals[22,22] = 0
Pvals = as.matrix(Pvals)

Pvals[Pvals>0.05]=0

network=graph_from_adjacency_matrix( Pvals, weighted=T, mode="undirected", diag=F)
plot(network)

# Example for how to change plotting parameters
#E(net)$color lets you access edge colors
#V(net)$color lets you access vertex colors

library(RColorBrewer)
xxx = brewer.pal(5,"Set2")
N1 = xxx[1] # Proteobacteria
N2 = xxx[2] # Actinobacteria
N3 = xxx[3] # Bacteroidetes
N4 = xxx[4] # Firmicutes
N5 = xxx[5] # Fusobacteria
# Old SigCorrs = c(0.352,0.341, 0.399,0.393,0.388,-0.383,-0.449,-0.390,0.437,0.419,0.384,-0.416,0.447)
SigCorrs = c(0.399, 0.393,0.388, -0.383, -0.449, -0.390, 0.437, 0.419,0.384,-0.416,0.447)
Colors = c("green","green","green","red","red","red", "green", "green", "green","red","green")
E(network)$color = Colors
V(network)$color = c(N1, N2, N3, N1, N4, N1, N4, N2, N5, N1, N4, N4, N1, N1, N1, N1, N1, N1, N1, N1, N3, N1)
plot(network, vertex.label.color = "black", vertex.size = 20, vertex.label.cex = 1.1, edge.curved = 0.1)

