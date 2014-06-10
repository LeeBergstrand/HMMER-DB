#!/usr/bin/env Rscript 
# Created by: Lee Bergstrand 
# Descript: Creates plot showing E-value distribution of HMMHits from HMMER-DB
#             
# Requirements: 
#-----------------------------------------------------------------------------
#Imports:
library(ggplot2)
library(scales)
library(RSQLite)

sqlite = dbDriver("SQLite")
HMMDB = dbConnect(sqlite, "/Users/lee/Data/SteriodHMMs/HMMDB.sqlite")

data = dbGetQuery(HMMDB, "SELECT HMM_Data.HMM_Family, count(HMM_Data.HMM_Family) as Count, substr(Organisms.Organism_Phylogeny,0,18) as Phylum
                          FROM HMM_Data, HMM_Hits, Organisms, Proteins
                          WHERE HMM_Data.HMM_Model = HMM_Hits.HMM_Model AND HMM_Hits.Protein_Accession = Proteins.Protein_Accession AND Proteins.Organism_Accession = Organisms.Organism_Accession
                          AND Organisms.Organism_Phylogeny LIKE 'Bacteria%'
                          GROUP BY HMM_Data.HMM_Family, Phylum
                          ORDER BY Phylum")

plotObj = ggplot(data, aes(x = HMM_Family, y = Count, colour = factor(Phylum), size = 3))
plotObj + geom_point() + facet_grid(Phylum ~ ., scales = "free")