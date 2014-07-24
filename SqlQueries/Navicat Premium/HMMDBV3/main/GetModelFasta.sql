SELECT DISTINCT Proteins.FASTA_Sequence
FROM HMM_Hits, Organisms, Proteins
WHERE HMM_Hits.Protein_Accession = Proteins.Protein_Accession
AND Proteins.Organism_Accession = Organisms.Organism_Accession
AND Organisms.Organism_Accession IN (
	/* SQL Inner Query 2: Selects organisms with an completeness (number of HMM families with hits) greater than 60% from organisms with hits for HsaC and KshA */
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
		AND HMM_Hits.HMM_Coverage >= 0.8 /* Select only hsaC hits with over 80% HMM coverage */
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
			AND HMM_Hits.HMM_Coverage >= 0.8 /* Select only KshA hits with over 80% HMM coverage */
		)
	)
	GROUP BY
		Organisms.Organism_Accession
	HAVING
		count(DISTINCT HMM_Data.HMM_Family) >= 13.8 /* Filters out organisms with an completeness (number of HMM families with hits) less than 60% */
	ORDER BY
		count(DISTINCT HMM_Data.HMM_Family) ASC
)
AND HMM_Hits.HMM_Model NOT LIKE '%HsaA%' /* Filters out Cyps */
AND HMM_Hits.HMM_Coverage >= 0.8
GROUP BY Proteins.Protein_Accession
HAVING max(HMM_Hits.HMM_Score)