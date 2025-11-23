### Building a phylogenetic tree
### Mathew REES
### Royal Botanic Garden Edinburgh
### mrees@rbge.org.uk
### 01/09/2025

### This script takes a list of species and builds a synthetic phylogenetic tree
### These trees are useful for large scale analyses but be wary when using them 
### for shallow taxonomic comparative analyses, as they will contain many polytomies.

library(tidyverse)
library(V.PhyloMaker) # Available here https://github.com/jinyizju/V.PhyloMaker
library(ape)
library(taxize) # version 0.10 and later required for updated API
library(stringdist)
library(tidystringdist)

# Here read in your species list
# This can be a single column with species, or ideally, it can be 3 columns 
# with species, genus and family, all lower case
# If you are workign with SEOSAW data, you should have these columns readily available
data <- read.csv("species_list.csv")

# The species names should be written without any special characters,
# like this "Brachystegia spiciformis"

# Make sure your list only has a single row per species, no duplicates
data <- data %>% 
  filter(!duplicated(species))

# Also make sure that all your names are spelled properly. No typos.
# Here we can check if there are any names that are very close in spelling
df <- tidy_comb_all(data$species)

df_dist <- tidy_stringdist(df)

df_dist %>% 
  filter(jw <= 0.1) # you can change which index you use as a filter and play with the threshold

# Make any changes accordingly

# if you don't have the genus or family columns, we can create these 
taxize_match <- gna_verifier(names = data$species, all_matches = F, data_sources = 196) # 196 here is the id for the World Flora Online 
# if you want to match names against a different database, use gna_data_sources()

# You can also check if your names are currently accepted or if they are synonyms
taxize_match %>% 
  dplyr::select(matchedCanonicalSimple, currentCanonicalSimple) %>% 
  mutate(dif = ifelse(matchedCanonicalSimple == currentCanonicalSimple, 1,0)) %>% 
  filter(dif == 0)

# We can now get the genus and family from NCBI
taxize_names <- tax_name(sci = gnr_names$matchedCanonicalFull, 
                         get = c("genus","family"), 
                         db="ncbi")

# Finalise your data object for the phylogenetic tree building function
data_phylo <- taxize_names %>% 
  dplyr::select(species = query, genus, family)

# Now we can run the V.Phylomaker function
# This uses the GBOTB.extended mega tree from Smith et al. 2018 and Magallon et al. 2014 by default
# This mega tree does contain some non-monophyletic genera

phylo <- phylo.maker(sp.list = data_phylo, scenarios = "S3")

# You can check how many of your species were already present in the tree (prune) and how many were grafted
head(phylo$species.list)
table(phylo$species.list$status)

# And check your phylogenetic tree
phylo$scenario.3
# Note here that the species names have received an underscore

# Write your phylogenetic tree to your local working folder
write.tree(phy = phylo$scenario.3, file = "my_phylogeny.tre")
