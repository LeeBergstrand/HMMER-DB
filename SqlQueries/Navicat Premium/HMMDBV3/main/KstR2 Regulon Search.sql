SELECT
	Proteins.Locus,
	Proteins.Protein_Accession,
	HMM_Hits.HMM_Model,
	HMM_Hits.HMM_Score,
	HMM_Hits.HMM_E_Value,
	HMM_Hits.HMM_Coverage
FROM
	HMM_Hits,
	Proteins
WHERE
	HMM_Hits.Protein_Accession = Proteins.Protein_Accession
AND Proteins.Locus IN (
	'Rv3548c',
	'Rv3549c',
	'Rv3550',
	'Rv3551',
	'Rv3552',
	'Rv3553',
	'Rv3556c',
	'Rv3557c',
	'Rv3559c',
	'Rv3560c',
	'Rv3561',
	'Rv3562',
	'Rv3563',
	'Rv3564',
	'Rv3565'
)
ORDER BY
	Proteins.Locus ASC,
	HMM_Hits.HMM_Score DESC