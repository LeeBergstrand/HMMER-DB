#!/usr/bin/env python 
# Created by: Lee Bergstrand 
# Descript: 
#
# Requirements: - This script requires the Biopython module: http://biopython.org/wiki/Download
#----------------------------------------------------------------------------------------

#===========================================================================================================
#Imports & Setup:
from Bio import SeqIO

#===========================================================================================================
# Functions:

# 1: When passed a sequence record object returns an array of fasta strings for each annotation.
def getProtienAnnotationFasta(seqRecord):
	fasta = []
	features = seqRecord.features # Each sequence has a list (called features) that stores seqFeature objects.
	for feature in features: # For each feature on the sequence
		if feature.type == "CDS": # CDS means coding sequence (These are the only feature we're interested in)
			featQualifers = feature.qualifiers # Each feature contains a dictionary called quailifiers which contains           
			                                   # data about the sequence feature (for example the translation)
			
			start  = int(feature.location.start) # Type-casting to int strips fuzzy < > characters.
			end    = int(feature.location.end)
			strand = feature.location.strand
			if strand == None:
				strand = "?"
			elif int(strand) < 0:
				strand = "-"
			elif int(strand) > 0:
				strand =  "+"
			else:
				strand = "?"
			
			location = "[" + str(start) + ":" + str(end) + "](" + strand + ")"
			
			# Gets the required qualifers. Uses featQualifers.get to return the quatifer or a default value if the quatifer
			# is not found. Calls strip to remove unwanted brackets and ' from quantifer before storing it as a string.
			protein_id = str(featQualifers.get('protein_id','no_protein_id')).strip('\'[]')
			if protein_id == 'no_protein_id':
				continue # Skips the iteration if protien has no id.
			protein_locus = str(featQualifers.get('locus_tag','no_locus_tag')).strip('\'[]')
			gene    = str(featQualifers.get('gene','no_gene_name')).strip('\'[]')
			product = str(featQualifers.get('product','no_product_name')).strip('\'[]')
			translated_protein = str(featQualifers.get('translation','no_translation')).strip('\'[]')
			fasta.append(">" + protein_id + " " + gene + "-" + product + " (Locus: " + protein_locus + ")" + " (Location: " + location + ")" + "\n" + translated_protein + "\n")
	FASTAString = "".join(fasta)
	return FASTAString
	
#------------------------------------------------------------------------------------------------------------
# 2: When passed a sequence record object returns a string containing the organisms info.
def getOrganismInfo(seqRecord):
	print ">> Extracting Organism Info..."
	OrganismID  = seqRecord.id
	description = seqRecord.description.replace(",", "")
	accessionType = ""
	if 'plasmid' in description.lower():
		accessionType = 'Plasmid'
	elif 'chromosome' in description.lower():
		accessionType = 'Chromosome'
	elif 'genome' in description.lower():
		accessionType = 'Chromosome'
	else:
		accessionType = 'Unknown'
	source = seqRecord.annotations['source'].replace(",", "")
	taxonomy = "_".join(seqRecord.annotations['taxonomy'])
	OrganismGenomeLength = len(seqRecord.seq)

	OrganismString = OrganismID + "," + accessionType + "," + description + "," + source + "," + taxonomy + "," + str(OrganismGenomeLength) + "\n"
	return OrganismString
#-------------------------------------------------------------------------------------------------------------