#!/usr/bin/env Rscript 
# Created by: Lee Bergstrand 
# Descript: Counts the number of HMM hits per HMM family per Phylum.
#             
# Requirements: ggplot2, scales, RSQLite
#-----------------------------------------------------------------------------
#Imports:
library(ggplot2)
library(scales)
library(RSQLite)

sqlite = dbDriver("SQLite")
HMMDB = dbConnect(sqlite, "/Users/lee/Data/SteriodHMMs/HMMDB.sqlite")

# Executes SQL query and loads results directly into a dataframe.
data = dbGetQuery(HMMDB, "/* SQL Query: Counts the number of distinct HMM families per organism source. */
                          SELECT
                              HMM_Data.HMM_Family,
	                            count(HMM_Data.HMM_Family) AS Count,
	                            substr(Organisms.Organism_Phylogeny,0,18) AS Phylum /* Feather the 18 to change the Phylogeny; Less:HigherTaxa <-> Greater:LowerTaxa */
                          FROM
	                            HMM_Data,
	                            HMM_Hits,
	                            Organisms,
	                            Proteins
                          WHERE
	                            HMM_Data.HMM_Model = HMM_Hits.HMM_Model
                          AND HMM_Hits.Protein_Accession = Proteins.Protein_Accession
                          AND Proteins.Organism_Accession = Organisms.Organism_Accession
                          AND Organisms.Organism_Phylogeny LIKE 'Bacteria%' /* Use this to filter by specific Phylogeny */
                          GROUP BY
	                            HMM_Data.HMM_Family,
	                            Phylum
                          ORDER BY
	                            Phylum")

# Plots Data as a bar graph.
plotObj = ggplot(data, aes(x = HMM_Family, y = Count, colour = factor(Phylum), size = 3))
plotObj + geom_point() + facet_grid(Phylum ~ ., scales = "free")