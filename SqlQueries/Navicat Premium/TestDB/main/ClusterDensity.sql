SELECT
	subQuery.HMM_Family,
	subQuery.Organism_Description,
	((subQuery.Protein_Relative_Start + subQuery.Protein_Relative_End) / 2) AS Protein_Relative_Center
FROM
	(
		SELECT DISTINCT
			HMM_Data.HMM_Family,
			Organisms.Organism_Description,
			(CAST(Proteins."End" AS float)/CAST(Organisms.Sequence_Length AS float) * 100) as Protein_Relative_End,
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
	) AS subQuery