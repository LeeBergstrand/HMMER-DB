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
	AND HMM_Data.HMM_Family IN ('hsaC', 'KshA')
)
AND Organisms.Organism_Phylogeny NOT LIKE '%Mycobacterium%'
AND Organisms.Organism_Phylogeny NOT LIKE '%Rhodococcus%'
GROUP BY
	Organisms.Source
HAVING
	count(DISTINCT HMM_Data.HMM_Family) >= 13.8
ORDER BY
	count(DISTINCT HMM_Data.HMM_Family) ASC