SELECT
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