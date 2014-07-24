SELECT
	count(subQuery.HMM_Family_Count) AS Completeness,
	count(subQuery.Source) AS Organism_Count,
	substr(subQuery.Phylogeny, 0, 20) AS Phylogeny
FROM
	(
		SELECT
			count(DISTINCT HMM_Data.HMM_Family) AS HMM_Family_Count,
			Organisms.Source AS Source,
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
		AND HMM_Data.HMM_Family NOT LIKE 'Fad%'
		AND HMM_Data.HMM_Family NOT LIKE 'Cyp%'
		AND HMM_Data.HMM_Family NOT LIKE 'EchA%'
		AND HMM_Data.HMM_Family NOT LIKE 'hsdA%'
		GROUP BY
			Organisms.Source
		ORDER BY
			HMM_Family_Count
	) AS subQuery
WHERE
	count(subQuery.HMM_Family_Count) >= 12
	Phylogeny LIKE 'Bacteria%'
GROUP BY
	Phylogeny
ORDER BY
	Phylogeny