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
HMMDB  = dbConnect(sqlite, "/Users/lee/Data/SteriodHMMs/OldDBs/HMMDBV4.sqlite") # Location of HMMER-DB Sqlite database.

# Executes SQL query and loads results directly into a dataframe. 
data = dbGetQuery(HMMDB, "/* SQL Outer Query: Wrapper for subquery one that calculates relative protein center from subquery one's relative protein start and stop. */
SELECT DISTINCT
  HMM_Data.HMM_Family,
                  Organisms.Organism_Description,
                  (( CAST (Proteins.'End' AS REAL) + CAST (Proteins.Start AS REAL)) / 2 ) AS Position
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
                  /* SQL Inner Query 1: Selects organisms from the subquery with a HMM_Coverage of at least 60% for key steriod degrading enzymes.*/
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
                  AND HMM_Data.HMM_Family IN (
                  'hsaA',
                  'hsaC',
                  'hsaD',
                  'hsaE',
                  'hsaF',
                  'hsaG',
                  'hsd4B',
                  'KshA',
                  'kstD'
                  )
                  AND HMM_Hits.HMM_Coverage >= 0.8
                  AND Organisms.Organism_Accession IN (
                  /* SQL Inner Query 2: Selects only organisms that hsaC hits with over 80% HMM coverage from organism with KshA*/
                  SELECT DISTINCT
                  Organisms.Organism_Accession
                  FROM
                  HMM_Hits,
                  Organisms,
                  Proteins
                  WHERE
                  HMM_Hits.Protein_Accession = Proteins.Protein_Accession
                  AND Organisms.Organism_Accession = Proteins.Organism_Accession
                  AND HMM_Hits.HMM_Model LIKE 'hsaC%' /* Selects organisms with a hit for HsaC*/
                  AND HMM_Hits.HMM_Coverage >= 0.8 /* Select only hsaC hits with over 80% HMM coverage */
                  AND Organisms.Organism_Accession IN (
                  /* SQL Inner Query 3: Selects only organisms that KshA hits with over 80% HMM coverage.*/
                  SELECT DISTINCT
                  Organisms.Organism_Accession
                  FROM
                  HMM_Hits,
                  Organisms,
                  Proteins
                  WHERE
                  HMM_Hits.Protein_Accession = Proteins.Protein_Accession
                  AND Organisms.Organism_Accession = Proteins.Organism_Accession
                  AND HMM_Hits.HMM_Model LIKE 'KshA%' /* Selects organisms with a hit for KshA*/
                  AND HMM_Hits.HMM_Coverage >= 0.8
                  AND Organisms.Organism_Description NOT LIKE 'Myco%' /* Selects non-Mycbacterium */
                  )
                  )
                  GROUP BY
                  Organisms.Organism_Description
                  HAVING
                  count(DISTINCT HMM_Data.HMM_Family) >= 5.4
                  ORDER BY
                  Organisms.Organism_Phylogeny ASC
                  )
                  AND HMM_Data.HMM_Family IN (
                  'hsaA',
                  'hsaC',
                  'hsaD',
                  'hsaE',
                  'hsaF',
                  'hsaG',
                  'hsd4B',
                  'KshA',
                  'kstD'
                  )
                  ORDER BY
                  Organisms.Organism_Phylogeny ASC,
                  Organisms.Organism_Description ASC,
                  Position ASC")



plotObj = ggplot(data, aes(x = Position, y = Organism_Description, color = factor(HMM_Family)))
plotObj + geom_point(alpha = 3/4) + # Slight alpha so one can visualize overlaping points better.
          ggtitle("Relative positions of proteins with HMM hits for organisms in the database.") + 
          xlab("Protein's postion on geneome (bp).") + ylab("Organism") + labs(colour = "Type of HMM hit") 
