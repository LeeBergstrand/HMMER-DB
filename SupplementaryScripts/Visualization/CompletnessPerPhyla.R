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

data = dbGetQuery(HMMDB, "SELECT avg(subQuery.HMM_Family_Count) as Completeness, count(subQuery.Source) as Organism_Count, substr(subQuery.Phylogeny,0,88) as Phylogeny
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
                          WHERE Phylogeny LIKE '%Actinobacteria%' AND subQuery.HMM_Family_Count >= 18
                          GROUP BY Phylogeny
                          ORDER BY Phylogeny")

print(data)  

plotObj = ggplot(data, aes(x = Phylogeny, y = Completeness, fill = Organism_Count))
plotObj + geom_bar(stat="identity") + coord_flip()

