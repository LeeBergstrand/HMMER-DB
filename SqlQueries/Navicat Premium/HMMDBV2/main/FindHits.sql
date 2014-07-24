SELECT DISTINCT HMM_Data.HMM_Family, 
	Proteins.Protein_Accession, 
	Proteins.Start, 
	Organisms.Organism_Accession, 
	Organisms.Sequence_Length, 
	Organisms.Organism_Description, 
	Organisms.Organism_Phylogeny
FROM HMM_Data, HMM_Hits, Organisms, Proteins
WHERE HMM_Data.HMM_Model = HMM_Hits.HMM_Model
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
ORDER BY Organisms.Organism_Description DESC, Proteins.Start ASC