SELECT DISTINCT HMM_Data.HMM_Family, 
	Proteins.Protein_Accession, 
	Proteins.Locus, 
	Proteins.Start, 
	Proteins."End", 
	Proteins.Strand, 
	Organisms.Source, 
	Organisms.Organism_Accession, 
	Organisms.Organism_Phylogeny, 
	Organisms.Accession_Type, 
	Organisms.Organism_Description
FROM HMM_Data, HMM_Hits, Organisms, Proteins
WHERE HMM_Data.HMM_Model = HMM_Hits.HMM_Model
AND HMM_Hits.Protein_Accession = Proteins.Protein_Accession
AND Proteins.Organism_Accession = Organisms.Organism_Accession
AND Organisms.Organism_Accession IN (
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
		AND HMM_Data.HMM_Family = 'hsaC'
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
			AND HMM_Data.HMM_Family IN ('KshA')
		)
	)
	/*AND Organisms.Organism_Phylogeny NOT LIKE '%Mycobacterium%'*/
	/*AND Organisms.Organism_Phylogeny NOT LIKE '%Rhodococcus%'*/
	GROUP BY
		Organisms.Organism_Accession
	HAVING
		count(DISTINCT HMM_Data.HMM_Family) >= 13.8
	ORDER BY
		count(DISTINCT HMM_Data.HMM_Family) ASC
)
AND HMM_Data.HMM_Family NOT LIKE '%FAD%'
AND HMM_Data.HMM_Family NOT LIKE '%EchA%'
AND HMM_Data.HMM_Family NOT LIKE '%hsd4A%'
AND HMM_Data.HMM_Family NOT LIKE '%Cyp%'
AND HMM_Hits.HMM_Coverage >= 0.8
ORDER BY Organisms.Accession_Type DESC, Organisms.Source ASC /*Organisms.Organism_Phylogeny DESC, Proteins.Locus ASC*/