#!/usr/bin/env Rscript 
# Created by: Lee Bergstrand 
# Descript: Creates a plot showing counts of hmm hits per HMM Family for each kingdom.
#             
# Requirements: ggplot2, scales, RSQLite
#-----------------------------------------------------------------------------
#Imports:
library(ggplot2)
library(scales)
library(RSQLite)

sqlite = dbDriver("SQLite")
HMMDB = dbConnect(sqlite, "/Users/lee/Data/SteriodHMMs/HMMDB.sqlite")

data = dbGetQuery(HMMDB, "SELECT
                              HMM_Data.HMM_Family,
	                            count(DISTINCT Organisms.Source) AS Organism_Count,
	                            substr(Organisms.Organism_Phylogeny,0,8) AS Kingdom
                          FROM
	                            HMM_Data,
	                            HMM_Hits,
	                            Organisms,
	                            Proteins
                          WHERE
	                            HMM_Data.HMM_Model = HMM_Hits.HMM_Model
                          AND HMM_Hits.Protein_Accession = Proteins.Protein_Accession
                          AND Proteins.Organism_Accession = Organisms.Organism_Accession
                          GROUP BY
	                            Kingdom,
	                            HMM_Data.HMM_Family
                          ORDER BY
	                            HMM_Data.HMM_Family,
	                            Organism_Count DESC")

# Plots Data as a bar graph.
plotObj = ggplot(data, aes(x = HMM_Family, y = Organism_Count, fill = Kingdom))
plotObj + geom_bar(stat = "identity", position = "identity") + coord_flip() + 
          theme(plot.background = element_rect(fill = "black"),
                title = element_text(colour = "white"),
                legend.title = element_text(colour = "black"),
                axis.title = element_blank(),
                axis.text = element_text(colour = "white"),
                axis.text.y = element_text(angle = -45, hjust = 1, size = rel(2)),
                axis.text.x = element_text(size = rel(2)),
                axis.ticks = element_line(colour = "white"))