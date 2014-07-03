#!/usr/bin/env Rscript 
# Created by: Lee Bergstrand 
# Descript: Creates plot showing position distribution of HMMHits from HMMER-DB
#             
# Requirements: ggplot2 and RSQLite
#-------------------------------------------------------------------------------

# Imports:
library(ggplot2)
library(RSQLite)

# Setting up database connection:
sqlite = dbDriver("SQLite")
HMMDB  = dbConnect(sqlite, "/Users/lee/Data/SteriodHMMs/HMMDBV2.sqlite") # Location of HMMER-DB Sqlite database.

# Executes SQL query and loads results directly into a dataframe. 
data = dbGetQuery(HMMDB, "/* SQL Outer Query: Wrapper for subquery one that calculates relative protein center from subquery one's relative protein start and stop.*/
                          SELECT
                            subQuery.HMM_Family,
                          	subQuery.Organism_Description,
                          	((subQuery.Protein_Relative_Start + subQuery.Protein_Relative_End) / 2) AS Protein_Relative_Center
                          FROM
                          	(
                              /* SQL Inner Query 1:  For each HMM_Hit, selects the HMM Family name and Organism with the hit and calculates the relative start and stop postions
                                                     of the hit's target protein. Does not select hits for FAD, EchA, hsd4A and Cyp series. The HMM coverage has to be greater
                                                     the 80% or the hit is filtered out. */
                          		SELECT DISTINCT
                          			HMM_Data.HMM_Family,
                          			Organisms.Organism_Description,
                          			(CAST(Proteins.'End' AS float)/CAST(Organisms.Sequence_Length AS float) * 100) as Protein_Relative_End,
                          			(CAST(Proteins.Start AS float)/CAST(Organisms.Sequence_Length AS float) * 100) as Protein_Relative_Start	
                          		FROM
                          			HMM_Data,
                          			HMM_Hits,
                          			Organisms,
                          			Proteins
                          		WHERE
                          			HMM_Data.HMM_Model = HMM_Hits.HMM_Model
                          		AND HMM_Hits.Protein_Accession = Proteins.Protein_Accession
                          		AND Proteins.Organism_Accession = Organisms.Organism_Accession
                          		AND Organisms.Organism_Accession IN (
                          			
                                    /* SQL Inner Query 2: Selects organisms with an completeness (number of HMM families with hits) greater than 60% from organisms with hits for HsaC and KshA*/
                                    SELECT
                          				Organisms.Organism_Accession
                          			FROM
                          				HMM_Data,
                          				HMM_Hits,
                          				Organisms,
                          				Proteins
                          			WHERE
                          				HMM_Data.HMM_Model = HMM_Hits.HMM_Model
                          			AND HMM_Hits.Protein_Accession = Proteins.Protein_Accession
                          			AND Proteins.Organism_Accession = Organisms.Organism_Accession
                          			AND Organisms.Organism_Accession IN (
                          				  
                                        /* SQL Inner Query 3: Selects organisms with hits for HsaC from organism with hits for KshA */
                                        SELECT DISTINCT
                          					Organisms.Organism_Accession
                          				FROM
                          					HMM_Data,
                          					HMM_Hits,
                          					Organisms,
                          					Proteins
                          				WHERE
                          					HMM_Data.HMM_Model = HMM_Hits.HMM_Model
                          				AND HMM_Hits.Protein_Accession = Proteins.Protein_Accession
                          				AND Organisms.Organism_Accession = Proteins.Organism_Accession
                          				AND HMM_Data.HMM_Family = 'hsaC' /* Selects organisms with a hit for HsaC*/
                          				AND Organisms.Organism_Accession IN (
                  
                            				/* SQL Inner Query 4: Selects organisms with hits for KshA */
                                            SELECT DISTINCT 
                          						Organisms.Organism_Accession
                          					FROM
                          						HMM_Data,
                          						HMM_Hits,
                          						Organisms,
                          						Proteins
                          					WHERE
                          						HMM_Data.HMM_Model = HMM_Hits.HMM_Model
                          					AND HMM_Hits.Protein_Accession = Proteins.Protein_Accession
                          					AND Organisms.Organism_Accession = Proteins.Organism_Accession
                          					AND HMM_Data.HMM_Family IN ('KshA') /* Selects organisms with a hit for KshA*/
                          				)
                          			)
                          			GROUP BY
                          				Organisms.Organism_Accession
                          			HAVING
                          				count(DISTINCT HMM_Data.HMM_Family) >= 13.8 /* Filters out organisms with an completeness (number of HMM families with hits) less than 60%*/
                          			ORDER BY
                          				count(DISTINCT HMM_Data.HMM_Family) ASC
                          		)
                          		AND HMM_Data.HMM_Family NOT LIKE '%FAD%' /* Filters out FADs */
                          		AND HMM_Data.HMM_Family NOT LIKE '%EchA%' /* Filters out EchAs */
                          		AND HMM_Data.HMM_Family NOT LIKE '%hsd4A%' /* Filters out hsd4As */
                          		AND HMM_Data.HMM_Family NOT LIKE '%Cyp%' /* Filters out Cyps */
                          		AND HMM_Hits.HMM_Coverage >= 0.8
                          	) AS subQuery")
print("Alah Spam")


plotObj = ggplot(data, aes(x = Protein_Relative_Center, y = Organism_Description, color = factor(HMM_Family)))
plotObj + geom_point(alpha = 3/4) + # Slight alpha so one can visualize overlaping points better.
          ggtitle("Relative positions of proteins with HMM hits for organisms in the database.") + 
          xlab("Protein's relative postion on geneome.") + ylab("Organism") + labs(colour = "Type of HMM hit") 
