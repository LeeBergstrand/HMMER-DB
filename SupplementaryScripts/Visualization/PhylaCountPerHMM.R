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
	                              substr(Organisms.Organism_Phylogeny,0,18) as Kingdom
                          FROM HMM_Data, HMM_Hits, Organisms, Proteins
                          WHERE HMM_Data.HMM_Model = HMM_Hits.HMM_Model 
                                AND HMM_Hits.Protein_Accession = Proteins.Protein_Accession 
                                AND Proteins.Organism_Accession = Organisms.Organism_Accession
                                AND Kingdom LIKE 'Bacteria%'
                                AND HMM_Data.HMM_Family NOT LIKE 'Fad%' 
                                AND HMM_Data.HMM_Family NOT LIKE 'Ech%'
                                AND HMM_Data.HMM_Family NOT LIKE 'hsd4A'
                                AND HMM_Data.HMM_Family NOT LIKE 'KshB'
                                AND HMM_Data.HMM_Family NOT LIKE 'Cyp142'
                                AND HMM_Data.HMM_Family NOT LIKE 'Cyp125'
                          GROUP BY Kingdom, HMM_Data.HMM_Family
                          ORDER BY HMM_Data.HMM_Family, Organism_Count DESC")

print(data)

plotObj = ggplot(data, aes(x = HMM_Family, y = Organism_Count, fill = Kingdom))
plotObj + geom_bar(stat = "identity", position = "identity") + coord_flip() + 
          theme(plot.background = element_rect(fill = "black"),
                title = element_text(colour = "white"),
                legend.title = element_text(colour = "black"),
                axis.title = element_text(colour = "white"),
                axis.text = element_text(colour = "white"),
                axis.text.y = element_text(angle = -45, hjust = 1),
                axis.ticks = element_line(colour = "white")) +
        xlab("HMM Family") +
        ylab("Organism Count") +
        ggtitle("Counts of organisms with HMM hits for specific proteins.")