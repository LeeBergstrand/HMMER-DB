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
from DBSetFunctions import getProtienAnnotationFasta
from DBSetFunctions import getOrganismInfo

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
#------------------------------------------------------------------------------------------------------------
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
