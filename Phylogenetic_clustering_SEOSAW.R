### Phylogenetic clustering
### Mathew REES
### Royal Botanic Garden Edinburgh
### mrees@rbge.org.uk
### 01/09/2025

### This script builds a phylogenetic distance matrix between pairs of plots,
### which can then be used for clustering
### It requires a phylogenetic tree and a site x species matrix (community matrix)
### Sites should be as rows and species as columns

library(tidyverse)
library(ape)
library(picante)
library(adiv)
library(ade4)
library(factoextra)
library(phyloregion)
library(cluster)
library(dendextend)

# First, read in your tree and your community matrix

tree <- read.tree(file = "my_tree.tree")
commat <- read.csv(file = "my_community_matrix")

# We want to make sure the tips of the phylogeny match the order of the species names
# in the community matrix
# Make sure that the format for writting the species names is the same in both objects,
# ie, underscores or other special characters

head(tree$tip.label)
head(colnames(commat))

com.phy <- match.phylo.comm(phy = tree, comm = commat)

## Now we can start building the phylogenetic dissimilarity matrix
## There are several ways we can do this:

## If our community matrix contains presence absence data, we can use the classic
## Sorensen dissimilarity index adapted to phylogenetic distance: Phylosor
## Or we can use the Jaccard dissimilarity index depending on whether we want to
## downweight the shared taxa. Look up the difference between these two indices to make your choice
## Good examples of this method are present in papers from Barnabas Daru and illustrated by the phyloregion package

## If your community matrix contains abundance data, we can use a Hellinger dissimilarity matrix.
## To do this we will first build a evolutionary PCA and then transform this PCA into a dissimilarity matrix
## This is the same approach used in Rees et al. 2025 Patterns and drivers of phylogenetic beta diversity

###############################################################################
## HELLINGER EVOPCA APPROACH ##
###############################################################################

## Assuming you are using SEOSAW data with abundance information, let's make the evoPCA
## if your commat included abundance data, leave abundance = TRUE
evocom <- evopcahellinger(phyl = com.phy$phy, comm = com.phy$comm, w = "evoab", scannf = F, nf = 20, abundance = T)

## you can visualise the PCA
fviz_pca(evocom)

## extract eigenvalues
get_eig(evocom)

# Now we can transform this Hellinger PCA into a Hellinger distance matrix
com.hel.dist <- dist.dudi(dudi = evocom)

###############################################################################
## PHYLOREGION APPROACH WITH SORENSEN INDEX
###############################################################################

betaphylosor <- phylobeta(x = dense2sparse(com.phy$comm), phy = com.phy$phy, index.family = "sorensen")

phylosor <- betaphylosor$phylo.beta.sor 
## note you can also focus on phylogenetic turnover if that is your interest
## see articles from Baselga 2010

################################################################################
## PHYLOGENETIC CLUSTERING
################################################################################

## You now have two dissimilarity matrices with which you can perform any type of clustering
## You should choose a clustering method appropriate to your question
## Here I will illustrate a simple hierarchical clustering with Ward's method

hc.hel <- hclust(d = com.hel.dist, method = "ward.D2")
plot(as.dendrogram(hc.hel))

hc.sor <- hclust(d = phylosor, method = "ward.D2")
plot(as.dendrogram(hc.sor))

## We now have to decide how many clusters there should be
## This is a topic beyond the scope of this script, but a few starters are provided below

optimal_k <- optimal_phyloregion(x = com.hel.dist, method = "ward.D2", k = 15)
paste0("otpimal number of clusters is ", optimal_k$optimal$k)

# Create a sequence of 1:15 to cut the tree at different heights
ags <- sapply(c(1:15), function(q) cutree(hc.hel, k=q)) ## make 1:15 clusters
ags <- data.frame(ags)
colnames(ags) <- paste("k", c(1:15), sep="")

## Each sample of your community matrix should now have a cluster assigned

# You can plot the optimal number of clusters, or any number of k
as.dendrogram(hc) %>%
  set("labels_col", value = palette.colors(n = optimal_k$optimal$k), k= optimal_k$optimal$k) %>%
  set("branches_k_color", value = palette.colors(n = optimal_k$optimal$k), k = optimal_k$optimal$k) %>%
  plot(horiz=TRUE, axes = F)
as.dendrogram(hc) %>% 
  rect.dendrogram(k = optimal_k$optimal$k, horiz = TRUE)
