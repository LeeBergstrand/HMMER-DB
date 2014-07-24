SELECT HMM_Data.HMM_Family, 
	count(HMM_Data.HMM_Family) as Count, 
	substr(Organisms.Organism_Phylogeny,0,19) as Phylum
FROM HMM_Data, HMM_Hits, Organisms, Proteins
WHERE HMM_Data.HMM_Model = HMM_Hits.HMM_Model 
	  AND HMM_Hits.Protein_Accession = Proteins.Protein_Accession 
	  AND Proteins.Organism_Accession = Organisms.Organism_Accession 
	  AND Organisms.Organism_Phylogeny LIKE "Bacteria_%"
GROUP BY HMM_Data.HMM_Family, Phylum 
ORDER BY Phylum