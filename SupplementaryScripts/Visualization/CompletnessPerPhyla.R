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
library(RSQLite.extfuns)

sqlite = dbDriver("SQLite")
HMMDB = dbConnect(sqlite, "/Users/lee/Data/SteriodHMMs/HMMDB.sqlite")
init_extensions(HMMDB)

data = dbGetQuery(HMMDB, "SELECT median(subQuery.HMM_Family_Count) as Completeness, count(subQuery.Source) as Organism_Count, substr(subQuery.Phylogeny,0,33) as Phylogeny
                          FROM
                          (  
                            SELECT count(DISTINCT HMM_Data.HMM_Family) as HMM_Family_Count, Organisms.Source as Source, Organisms.Organism_Phylogeny as Phylogeny
                          	FROM HMM_Data, HMM_Hits, Organisms, Proteins
                          	WHERE HMM_Data.HMM_Model = HMM_Hits.HMM_Model 
                                  AND HMM_Hits.Protein_Accession = Proteins.Protein_Accession 
                                  AND Proteins.Organism_Accession = Organisms.Organism_Accession
                          	GROUP BY Organisms.Source
                          	ORDER BY HMM_Family_Count
                          ) as subQuery
                          WHERE Phylogeny LIKE 'Bacteria%'
                          GROUP BY Phylogeny
                          HAVING Completeness >= 11.5
                          ORDER BY Phylogeny")

print(data)  

plotObj = ggplot(data, aes(x = Phylogeny, y = Completeness))
plotObj + geom_bar(stat="identity") + coord_flip() +
          theme(plot.background = element_rect(fill = "black"),
                title = element_text(colour = "white"),
                legend.title = element_text(colour = "black"),
                legend.position = "none",
                axis.title = element_blank(),
                axis.text = element_text(colour = "white"),
                axis.text.y = element_text(size = rel(2)),
                axis.text.x = element_text(size = rel(2)),
                axis.ticks = element_line(colour = "white"))
