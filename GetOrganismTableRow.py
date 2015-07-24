#!/usr/bin/env python 
# Created by: Lee Bergstrand 
# Descript: Appends a csv file with the organism info.
#
# Requirements: - This script requires the Biopython module: http://biopython.org/wiki/Download
#
# Usage: GenbankToFASTAandOrganismTable.py <organism.gbk> 
# Example: GenbankToFASTAandOrganismTable.py ecoli.gbk
# ----------------------------------------------------------------------------------------
# ===========================================================================================================
# Imports & Setup:

import sys
from Bio import SeqIO
from os import path


# ===========================================================================================================
# Functions:
# ------------------------------------------------------------------------------------------------------------
# 1: Checks if in proper number of arguments are passed gives instructions on proper use.
def argsCheck(numArgs):
	if len(sys.argv) < numArgs or len(sys.argv) > numArgs:
		print("Extracts annotation information and organism info from a genbank file.")
		print("By Lee Bergstrand\n")
		print("Usage: " + sys.argv[0] + " <organism.gbk>")
		print("Examples: " + sys.argv[0] + " ecoli.gbk")
		sys.exit(1)  # Aborts program. (exit(1) indicates that an error occurred)


# ------------------------------------------------------------------------------------------------------------
# 2: When passed a sequence record object returns a string containing the organisms info.
def getOrganismInfo(seqRecord):
	print(">> Extracting Organism Info...")
	OrganismID = seqRecord.id
	description = seqRecord.description.replace(",", "")
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

	OrganismString = OrganismID + "," + accessionType + "," + description + "," + source + "," + taxonomy + "," + str(
		OrganismGenomeLength) + "\n"
	return OrganismString


# -------------------------------------------------------------------------------------------------------------

# Main program code:

# House keeping...
argsCheck(2)  # Checks if the number of arguments are correct.

# Stores file one for input checking.
print(">> Starting up...")
OrganismFile = sys.argv[1]
OrganismID = path.split(OrganismFile)[1].rstrip(".gbk")

# File extension checks
print(">> Performing file extension checks...")
if not OrganismFile.endswith(".gbk"):
	print("[Warning] " + OrganismFile + " may not be a genbank file!")

# Read in genbank file as a sequence record object.
try:
	print(">> Opening " + OrganismFile)
	handle = open(OrganismFile, "rU")
	try:
		record = SeqIO.read(handle, "genbank")
	except ValueError as error:
		print("Error has occured while parsing " + OrganismFile + "!")
		print(error)
		sys.exit(1)
	handle.close()
except IOError:
	print("Failed to open " + OrganismFile)
	sys.exit(1)

print(">> Extracting Protein Annotations...")
OrganismString = getOrganismInfo(record)  # Creates a comma separated string with organism info.

# Appends to organism CSV.	
try:
	print(">> Writing to organism info to CSV file...")
	FASTAWriter = open("OrganismDB.csv", "w")
	FASTAWriter.write(OrganismString)
	FASTAWriter.close()
except IOError:
	print("Failed to open OrganismDB.csv")
	sys.exit(1)
print(">> Done")
