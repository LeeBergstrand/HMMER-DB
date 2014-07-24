SELECT
	count(subQuery.Source) AS Organism_Count,
	subQuery.HMM_Family_Count AS Proteins_Count,
	substr(subQuery.Phylogeny, 0, 8) AS Kingdom
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
		GROUP BY
			Organisms.Source
		ORDER BY
			HMM_Family_Count
	) AS subQuery
GROUP BY
	Kingdom,
	Proteins_Count