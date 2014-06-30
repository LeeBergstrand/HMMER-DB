#!/usr/bin/env Rscript 
# Created by: Lee Bergstrand 
# Descript: Creates a plot that shows the completeness (number of hmm familes with hits) per phylum.
#             
# Requirements: ggplot2, scales, RSQLite, RSQLite.extfuns
#-----------------------------------------------------------------------------
#Imports:
library(ggplot2)
library(scales)
library(RSQLite)
library(RSQLite.extfuns)

sqlite = dbDriver("SQLite")
HMMDB = dbConnect(sqlite, "/Users/lee/Data/SteriodHMMs/HMMDB.sqlite")
init_extensions(HMMDB)

# Executes SQL query and loads results directly into a dataframe.
# SQL Inner Query:  Counts the number of distinct HMM families per organism source.
# SQL Outers Query: Takes the median of counts above counts per user definded Phylogeny substring.
data = dbGetQuery(HMMDB, "SELECT
                              median (subQuery.HMM_Family_Count) AS Completeness,
                              substr(subQuery.Phylogeny, 0, 300) AS Phylogeny, /* Feather the 33 to change the Phylogeny; Less:HigherTaxa <-> Greater:LowerTaxa */
                              substr(subQuery.Phylogeny, 25, 15) AS Phylum /* Feather the 25/15 to change the what taxa are coloured and on the legend */
                          FROM
                            (
                              SELECT
                                  count(DISTINCT HMM_Data.HMM_Family) AS HMM_Family_Count,
                                  Organisms.Organism_Phylogeny AS Phylogeny
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
                                  Organisms.Source
                              ORDER BY
                                  HMM_Family_Count
                            ) AS subQuery
                          WHERE
                              Phylogeny LIKE '%Actinobacteria%' /* Use this to filter by specific Phylogeny */
                          GROUP BY
                              Phylogeny
                          HAVING
                              Completeness >= 13.8 /* Use this to filter by specific completeness (Depends on dataset) */
                          ORDER BY
                              Phylogeny")

# Plots Data as a bar graph.
plotObj = ggplot(data, aes(x = Phylogeny, y = Completeness, fill = Phylum))
plotObj + geom_bar(stat="identity") + coord_flip() +
          theme(plot.background = element_rect(fill = "black"),
                title = element_text(colour = "white"),
                legend.title = element_text(colour = "black"),
                axis.title = element_blank(),
                axis.text = element_text(colour = "white"),
                axis.text.y = element_text(size = rel(1.5)),
                axis.text.x = element_text(size = rel(1.5)),
                axis.ticks = element_line(colour = "white"))