#!/usr/bin/env python 
# Created by: Lee Bergstrand 
# Descript: A program that extracts the protiens annotations from a genbank file and as well as some 
#			information about the organism in the file. Stores the protein annotations as a Fasta. 
#			Appends a csv file with the organism info.
#
# Requirements: - This script requires the Biopython module: http://biopython.org/wiki/Download
#
# Usage: GenbankToFASTAandOrganismTable.py <organism.gbk> 
# Example: GenbankToFASTAandOrganismTable.py ecoli.gbk
#----------------------------------------------------------------------------------------
#===========================================================================================================
#Imports & Setup:
	
import sys
from Bio import SeqIO
from os import path

#===========================================================================================================
# Functions:

# 1: Checks if in proper number of arguments are passed gives instructions on proper use.
def argsCheck(numArgs):
	if len(sys.argv) < numArgs or len(sys.argv) > numArgs:
		print "Extracts annotation information and organism info from a genbank file."
		print "By Lee Bergstrand\n"
		print "Usage: " + sys.argv[0] + " <organism.gbk>"
		print "Examples: " + sys.argv[0] + " ecoli.gbk"
		sys.exit(1) # Aborts program. (exit(1) indicates that an error occurred)
		
#-------------------------------------------------------------------------------------------------------------
# 2: When passed a sequence record object returns an array of fasta strings for each annotation.
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
# 3: When passed a sequence record object returns a string containing the organisms info.
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

# Main program code:

# House keeping...
argsCheck(2) # Checks if the number of arguments are correct.

# Stores file one for input checking.
print ">> Starting up..."
OrganismFile = sys.argv[1]
OrganismID = path.split(OrganismFile)[1].rstrip(".gbk")
FastaFile = OrganismID + ".faa"

# File extension checks
print ">> Performing file extention checks..."
if not OrganismFile.endswith(".gbk"):
	print "[Warning] " + OrganismFile + " may not be a genbank file!"

# Read in genbank file as a sequence record object.
try:
	print ">> Opening " + OrganismFile
	handle = open(OrganismFile, "rU")
	try:
		record = SeqIO.read(handle, "genbank")
	except ValueError as error:
		print "Error has occured while parsing " + OrganismFile + "!"
		print error
		print error.traceback
		sys.exit(1)
	handle.close()
except IOError:
	print "Failed to open " + OrganismFile
	sys.exit(1)
	
print ">> Extracting Protein Annotations..."
FASTA = getProtienAnnotationFasta(record) # Creates a string containing all protein annotations in the gbk file.
OrganismString = getOrganismInfo(record) # Creates a comma seperated string with organism info.

# Writes annotations to FASTA file.
try:
	print ">> Writing FASTA File..."
	FASTAWriter = open(FastaFile, "w")
	FASTAWriter.write(FASTA)
	FASTAWriter.close()
except IOError:
	print "Failed to open " + FastaFile
	sys.exit(1)

# Appends to organism CSV.	
try:
	print ">> Writing to organism info to CSV file..."
	FASTAWriter = open("OrganismDB.csv", "a")
	FASTAWriter.write(OrganismString)
	FASTAWriter.close()
except IOError:
	print "Failed to open " + CSVFile
	sys.exit(1)
print ">> Done"
