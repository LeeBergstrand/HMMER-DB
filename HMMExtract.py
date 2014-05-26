#!/usr/bin/env python 
# Created by: Lee Bergstrand 
# Descript: A program that extracts the protiens annotations from a genbank file and searches these 
#           annotations using HMMsearch and an HMM file. It then stores hits in sqlite database. 
#
# Requirements: - This script requires the Biopython module: http://biopython.org/wiki/Download
#               - This script requires HMMER 3.0 or later.
#               - This script requires sqlLite3 or later.
#
# Usage: HMMExtract.py <organism.gbk> <hmm.hmm> <sqldb.sqlite>
# Example: HMMExtract.py ecoli.gbk helicase.hmm helicasedb.sqlite
#----------------------------------------------------------------------------------------
#===========================================================================================================
#Imports & Setup:
	
import sys
import subprocess
import sqlite3
import cStringIO
import re
from os import path
from Bio import SeqIO 
from multiprocessing import cpu_count

import time

processors = cpu_count() # Gets number of processor cores for HMMER.

#===========================================================================================================
# Functions:

# 1: Checks if in proper number of arguments are passed gives instructions on proper use.
def argsCheck(numArgs):
	if len(sys.argv) < numArgs or len(sys.argv) > numArgs:
		print "Sequence Downloader"
		print "By Lee Bergstrand\n"
		print "Usage: " + sys.argv[0] + " <organism.gbk> <hmm.hmm> <sqldb.sqlite>"
		print "Examples: " + sys.argv[0] + " ecoli.gbk helicase.hmm helicasedb.sqlite"
		sys.exit(1) # Aborts program. (exit(1) indicates that an error occurred)
#------------------------------------------------------------------------------------------------------------
# 2: Runs HMMER with settings specific for extracting subject sequences.
def runHMMSearch(FASTA, HMMERDBFile):
	Found16S = True
	process = subprocess.Popen(["hmmsearch", "--acc", "--cpu", str(processors), HMMERDBFile, "-"], stdin = subprocess.PIPE, stdout = subprocess.PIPE, bufsize = 1)
	stdout, error = process.communicate(FASTA) # This returns a list with both stderr and stdout. Only want stdout which is first element.
	if error:
		print error
		sys.exit(1)
	else:
		return stdout
#------------------------------------------------------------------------------------------------------------
# 3: When passed a sequence record object returns an array of fasta strings for each annotation.
def getProtienAnnotationFasta(seqRecord):
	fasta = {}
	features = seqRecord.features # Each sequence has a list (called features) that stores seqFeature objects.
	for feature in features: # For each feature on the sequence
		if feature.type == "CDS": # CDS means coding sequence (These are the only feature we're interested in)
			featQualifers = feature.qualifiers # Each feature contains a dictionary called quailifiers which contains           
			                                   # data about the sequence feature (for example the translation)
			
			# Gets the required qualifers. Uses featQualifers.get to return the quatifer or a default value if the quatifer			# is not found. Calls strip to remove unwanted brackets and ' from quantifer before storing it as a string.
			protein_id = str(featQualifers.get('protein_id','no_protein_id')).strip('\'[]')
			if protein_id == 'no_protein_id':
				continue # Skips the iteration if protien has no id.
			gene    = str(featQualifers.get('gene','no_gene_name')).strip('\'[]')
			product = str(featQualifers.get('product','no_product_name')).strip('\'[]')
			translated_protein = str(featQualifers.get('translation','no_translation')).strip('\'[]')
			fasta[protein_id] = (">" + protein_id + " " + gene + "-" + product + "\n" + translated_protein + "\n")
	return fasta
#------------------------------------------------------------------------------------------------------------
# 4: Parses HMM searches text output and generates a two dementional array of the domain alignments results.
def parseHmmsearchResults(HMMResults, HMMName, HMMLength):
	HitRowRex = re.compile("^\s*\d\s*((\?)|(\!))\s*")
	HMMResults = HMMResults.split(">>") # Splits output at domain alignments.
	del HMMResults[0] # Deletes stuff at top of text output which would be the first element after the split.
	HMMResults = [x.split("Alignments")[0] for x in HMMResults] # Removes detailed per alignment info.
	HMMResultsCleaned = []
	for proteinResult in HMMResults:
		proteinResult = proteinResult.splitlines()
		TargetProtein = proteinResult[0].split()[0] # Records protein accession from first line
		for row in proteinResult:
			if HitRowRex.match(row): # If row is a domain table line.
				row = row.split()
				score   = float(row[2])
				evalue  = float(row[5])
				hmmfrom = int(row[6])
				hmmto   = int(row[7])
				alifrom = int(row[9])
				alito   = int(row[10])
				hmmCoverage = ((float((hmmto - hmmfrom)))/float((HMMLength)))
				DomainRow = [TargetProtein, HMMName, score, evalue, hmmfrom, hmmto, alifrom, alito, hmmCoverage]
				HMMResultsCleaned.append(DomainRow)
	return HMMResultsCleaned
#------------------------------------------------------------------------------------------------------------
# 6: Fitres HMM Hits.
def fitreHMMHitTable(HMMHitTable):
	#HMMHitTable[:] = [row for row in HMMHitTable if row[3] < float("1e-30")] # Filtres by E-value.
	
	FiltredHMMHitTable = list(HMMHitTable)
	for i in range(0,(len(HMMHitTable)-1)):
		RowOne = FiltredHMMHitTable[i]    # Current Row in hit table.
		RowTwo = FiltredHMMHitTable[i+1]  # Row below.
		if(RowOne[0] == RowTwo[0]): # If they have the same targe protein.
			AlignmentOneLength = RowOne[-2] - RowOne[-3] # RowOne AliTo - AliFrom
			AlignmentTwoLength = RowTwo[-2] - RowTwo[-3] # RowTwo AliTo - AliFrom
			Overlap = RowOne[-2] - RowTwo[-3]  # RowOne AliTo -  RowTwo AliFrom
			print Overlap
			if (Overlap > 0): # If there is overlap
				print "Overlap"
				# If the overlap is greater than 50% of either alignment.
				if((((float(Overlap)/float(AlignmentOneLength)) > 0.5) or ((float(Overlap)/float(AlignmentOneLength)) > 0.5))): 
				print "Is greater than 50"
						HMMHitTable.remove(RowOne) # The alignment with the lower evalue will allways be in the second row since 
												    # HMMSearch always formats the output with multiple domains highest to lowest E-value.
				else:
					print "Is greater than 50"
			else:
				print "No Overlap"
	
	#HMMHitTable[:] = [row for row in HMMHitTable if row[-1] > 0.3] # Filtres by Query Coverage.
	return	HMMHitTable		
#------------------------------------------------------------------------------------------------------------
# 6: Creates list of hits protien FASTAs.
def getHitProteins(HMMHitTable, AnnotationFASTADict):
	HitProteins = []
	for row in HMMHitTable:
		HitProteins.append(AnnotationFASTADict[row[0]])
	return HitProteins
#===========================================================================================================
# Main program code:
	
# House keeping...
argsCheck(4) # Checks if the number of arguments are correct.

# Stores file one for input checking.
print ">> Starting up..."
OrganismFile = sys.argv[1]
HMMFile      = sys.argv[2]
sqlFile      = sys.argv[3]

HMMName = path.split(HMMFile)[1].rstrip(".hmm")

# File extension checks
print ">> Performing file extention checks..."
if not OrganismFile.endswith(".gbk"):
	print "[Warning] " + OrganismFile + " may not be a genbank file!"
if not HMMFile.endswith(".hmm"):
	print "[Warning] " + HMMFile + " may not be a HMM file!"
if not sqlFile.endswith(".sqlite"):
	print "[Warning] " + sqlFile + " may not be a sqlite file!"

# Read in genbank file as a sequence record object.
try:
	print ">> Opening Genbank File..."
	handle = open(OrganismFile, "rU")
	record = SeqIO.read(handle, "genbank")
	handle.close()
except IOError:
	print "Failed to open " + OrganismFile
	sys.exit(1)

# Gets HMM Length.
try:
	HMMLengthRegex = re.compile("^LENG\s*\d*$")
	print ">> Opening HMM File..."
	with open(HMMFile, "rU") as inFile:
		currentLine = inFile.readline()
		while not HMMLengthRegex.match(currentLine):
			currentLine = inFile.readline()
		
		HMMLength = int(currentLine.split()[1])
except IOError:
	print "Failed to open " + HMMFile
	sys.exit(1)
	
print ">> Extracting Protein Annotations..."
AnnotationFASTADict = getProtienAnnotationFasta(record) # Creates a dictionary containing all protein annotations in the gbk file.
print ">> Extracting Organism Info..."

OrganismInfo = [record.annotations['source'], record.annotations['taxonomy'][0], record.annotations['taxonomy']]
FASTAString = "".join(AnnotationFASTADict.values()) # Saves these annotations to a string.

#HMMResults = runHMMSearch(FASTAString, HMMFile) # Runs hmmsearch.

#HMMHitTable = parseHmmsearchResults(HMMResults, HMMName, HMMLength) # Parses hmmsearch results into a two dimensional array.

HMMHitTable = [['YP_705149.1', 'Cyp125(GramPos)', 1.3, 42.0, 99, 132, 104, 137, 0.08088235294117647]]
HMMHitTable.append(['YP_705149.1', 'Cyp125(GramPos)', 30.4, 5.9e-08, 214, 376, 230, 420, 0.39705882352941174])

HMMHitTable.append(['YP_703037.1', 'Cyp125(GramPos)', 2.7, 16.0, 100, 130, 107, 236, 0.07352941176470588])
HMMHitTable.append(['YP_703037.1', 'Cyp125(GramPos)', 11.5, 0.035, 213, 275, 231, 293, 0.15196078431372548])
HMMHitTable.append(['YP_703037.1', 'Cyp125(GramPos)', 18.1, 0.00034, 286, 376, 323, 422, 0.22058823529411764])

HMMHitTable.append(['YP_704615.1', 'Cyp125(GramPos)', 6.1, 1.4, 198, 273, 196, 278, 0.18382352941176472])
HMMHitTable.append(['YP_704615.1', 'Cyp125(GramPos)', 26.7, 7.9e-07, 283, 377, 308, 411, 0.23039215686274508])

HMMHitTable.append(['YP_704571.1', 'Cyp125(GramPos)', -3.4, 1100.0, 104, 148, 136, 180, 0.10784313725490197])
HMMHitTable.append(['YP_704571.1', 'Cyp125(GramPos)', 14.5, 0.0041, 213, 377, 270, 464, 0.4019607843137255])



for i in HMMHitTable:
	print i 
HMMHitTable = fitreHMMHitTable(HMMHitTable)

# HitProtienFASTAs = getHitProteins(HMMHitTable, AnnotationFASTADict) # Gets hit protein FASTAs.