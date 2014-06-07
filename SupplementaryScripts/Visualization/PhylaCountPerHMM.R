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

data = dbGetQuery(HMMDB, "SELECT HMM_Data.HMM_Family, 
                                count(DISTINCT Organisms.Source) as Organism_Count, 
	                              substr(Organisms.Organism_Phylogeny,0,8) as Kingdom
                          FROM HMM_Data, HMM_Hits, Organisms, Proteins
                          WHERE HMM_Data.HMM_Model = HMM_Hits.HMM_Model 
                                AND HMM_Hits.Protein_Accession = Proteins.Protein_Accession 
                                AND Proteins.Organism_Accession = Organisms.Organism_Accession
                          GROUP BY Kingdom, HMM_Data.HMM_Family
                          ORDER BY HMM_Data.HMM_Family, Organism_Count DESC")

print(data)  

plotObj = ggplot(data, aes(x = HMM_Family, y = Organism_Count, fill = Kingdom))
plotObj + geom_bar(stat = "identity", position = "identity") + coord_flip()