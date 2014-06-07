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

data = dbGetQuery(HMMDB, "SELECT count(subQuery.Source) as Organism_Count, subQuery.HMM_Family_Count as Proteins_Count, substr(subQuery.Phylogeny,0,18) as Phyla
                          FROM
                          (  
                          	SELECT count(DISTINCT HMM_Data.HMM_Family) as HMM_Family_Count, Organisms.Source as Source, Organisms.Organism_Phylogeny as Phylogeny
                          	FROM HMM_Data, HMM_Hits, Organisms, Proteins
                          	WHERE HMM_Data.HMM_Model = HMM_Hits.HMM_Model AND HMM_Hits.Protein_Accession = Proteins.Protein_Accession AND Proteins.Organism_Accession = Organisms.Organism_Accession
                          	GROUP BY Organisms.Source
                          	ORDER BY HMM_Family_Count
                          ) as subQuery
                          WHERE Phyla LIKE 'Bacteria%'
                          GROUP BY Phyla, Proteins_Count
                          ORDER BY Proteins_Count, Organism_Count DESC")

print(data)  

plotObj = ggplot(data, aes(x = Proteins_Count, y = Organism_Count, fill = Phyla))
plotObj + geom_bar(stat = "identity", position = "dodge") + labs(title = 'Number of organisms with a specific number distince HMM hits per Phyla.') +
          scale_x_continuous(breaks = round(seq(min(data$Proteins_Count), max(data$Proteins_Count), by = 1))) +  facet_wrap(~ Phyla, scales = "free_x", nrow = 8, ncol = 4) +
          scale_y_continuous(breaks = round(seq(0, max(data$Organism_Count), by = 25))) + theme_bw() + theme(legend.position = "none") + 
          ylab("Number of organisms with a specific completeness") +
          xlab("Completeness (Number of distinct HMM Hits)")