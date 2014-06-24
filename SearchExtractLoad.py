#!/usr/bin/env python 
# Created by: Lee Bergstrand 
# Descript: A program that extracts the protein annotations from a fasta file and searches these 
#           annotations using HMMsearch and an HMM file. It then stores hits along with organism
#           information (gathered from a csv file) in a sqlite3 database. 
#
# Requirements: - This script requires the Biopython module: http://biopython.org/wiki/Download
#               - This script requires HMMER 3.1 or later.
#               - This script requires sqlLite3 or later.
#
# Usage: HMMExtract.py <organism.faa> <organisms.csv> <hmm.hmm> <sqldb.sqlite>
# Example: HMMExtract.py ecoli.faa bacteria.csv helicase.hmm helicasedb.sqlite
#----------------------------------------------------------------------------------------
#===========================================================================================================
#Imports & Setup:
	
import sys
import subprocess
import sqlite3
import csv
import re
import hashlib
from os import path
from Bio import SeqIO 
from multiprocessing import cpu_count

import time # Dev. Import

processors = cpu_count() # Gets number of processor cores for HMMER.

# Regex's
LocusRegex = re.compile("\(Locus:\s\S*\)")
LocationRegex = re.compile("\(Location:\s\[(\S*)\:(\S*)\]\((\S)\)\)")
#===========================================================================================================
# Functions:

# 1: Checks if in proper number of arguments are passed gives instructions on proper use.
def argsCheck(numArgs):
	if len(sys.argv) < numArgs or len(sys.argv) > numArgs:
		print "Hmmer hit cataloguer."
		print "By Lee Bergstrand\n"
		print "Usage: " + sys.argv[0] + "  <organism.faa> <organisms.csv> <hmm.hmm> <sqldb.sqlite>"
		print "Examples: " + sys.argv[0] + "  <organism.faa> <organisms.csv> <hmm.hmm> <sqldb.sqlite>"
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
# 3: Creates a python dictionary (hash table) that contains the the fasta for each protien in the proteome.
def createProteomeHash(ProteomeFile):
	ProteomeHash = dict() 
	try:
		handle = open(ProteomeFile, "rU")
		proteome = SeqIO.parse(handle, "fasta")
		for record in proteome:
			ProteomeHash.update({ record.id : record.format("fasta") })
		handle.close()
	except IOError:
		print "Failed to open " + ProteomeFile
		exit(1)
	return ProteomeHash
#------------------------------------------------------------------------------------------------------------
# 4: Parses HMM searches text output and generates a two dementional array of the domain alignments results.
def parseHmmsearchResults(HMMResults, HMMName, HMMLength):
	HitRowRegex = re.compile("^\s*\d\s*((\?)|(\!))\s*")
	HMMResults = HMMResults.split(">>") # Splits output at domain alignments.
	del HMMResults[0] # Deletes stuff at top of text output which would be the first element after the split.
	HMMResults = [x.split("Alignments")[0] for x in HMMResults] # Removes detailed per alignment info.
	HMMResultsCleaned = []
	for proteinResult in HMMResults:
		proteinResult = proteinResult.splitlines()
		TargetProtein = proteinResult[0].split()[0] # Records protein accession from first line
		for row in proteinResult:
			if HitRowRegex.match(row): # If row is a domain table line.
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
def filterHMMHitTable(HMMHitTable):

	HMMHitTable = [row for row in HMMHitTable if row[3] < float("1e-30")] # Filtres by E-value.
	
	i = 0
	while i < (len(HMMHitTable)-1):
		RowOne = HMMHitTable[i]    # Current Row in hit table.
		RowTwo = HMMHitTable[i+1]  # Row below.
		if(RowOne[0] == RowTwo[0]): # If they have the same targe protein.
			AlignmentOneLength = RowOne[-2] - RowOne[-3] # RowOne AliTo - AliFrom
			AlignmentTwoLength = RowTwo[-2] - RowTwo[-3] # RowTwo AliTo - AliFrom
			Overlap = RowOne[-2] - RowTwo[-3]  # RowOne AliTo -  RowTwo AliFrom
			if (Overlap > 0): # If there is overlap...
				# If the overlap is greater than 50% of either alignment.
				if((((float(Overlap)/float(AlignmentOneLength)) > 0.5) or ((float(Overlap)/float(AlignmentTwoLength)) > 0.5))):
					if RowOne[3] < RowTwo[3]: # If row one has a lower E-value than row two remove row two else remove row one.
						HMMHitTable.remove(RowTwo)
					else:
						HMMHitTable.remove(RowOne)
					i = i - 1 # Resets list index.
		i += 1
			
	HMMHitTable = [row for row in HMMHitTable if row[-1] > 0.3] # Filtres by Query Coverage.
	return	HMMHitTable		
#------------------------------------------------------------------------------------------------------------
# 6: Creates list of hits protien FASTAs.
def getHitProteins(HMMHitTable, AnnotationFASTADict, OrganismName):
	HitProteins = []
	for row in HMMHitTable:
		ProteinAccession = row[0]
		ProteinFASTA = AnnotationFASTADict[ProteinAccession]

		Locus = LocusRegex.search(ProteinFASTA).group(0)
		Locus = Locus.split()[1].rstrip(")")
		try:
			LocationData = LocationRegex.search(ProteinFASTA)
			Start  = LocationData.group(1)
			End    = LocationData.group(2)
			Strand = LocationData.group(3)
			ProteinData = [ProteinAccession, OrganismName, Locus, Start, End, Strand, ProteinFASTA]
			HitProteins.append(ProteinData)
		except AttributeError as Error:
			print row
			print ProteinFASTA
			print LocationData
			print "This is the organism: ", OrganismName
			print "The AttributeError was ", str(Error)
			sys.exit(1)
	return HitProteins
#-----------------------------------------------------------------------------------------------------------
# 7: Inserts organism info into DB.
def insertOrganismInfo(cursor, OrganismInfo):
	cursor.execute('''INSERT OR REPLACE INTO Organisms(Organism_Accession, Accession_Type, Organism_Description, Source, Organism_Phylogeny, Sequence_Length)
			          VALUES(?,?,?,?,?,?)''', OrganismInfo)
#-----------------------------------------------------------------------------------------------------------
# 8: Inserts protein info into DB.
def insertProteins(cursor, HitProteins):
	for protein in HitProteins:
		cursor.execute('''INSERT OR REPLACE INTO Proteins(Protein_Accession,Organism_Accession,Locus,Start,End,Strand,FASTA_Sequence)
		                  VALUES(?,?,?,?,?,?,?)''', protein)
#-----------------------------------------------------------------------------------------------------------
# 9: Inserts hits into DB and creates md5 hash for primary key.
def insertHits(cursor, HMMHitTable):
	HitHash = hashlib.md5()
	for hit in HMMHitTable:
		HitHash.update("".join([str(i) for i in HMMHitTable])) # Converts every list element to a string, joins them and updates Hash.
		hit = [HitHash.hexdigest()] + hit
		cursor.execute('''INSERT OR REPLACE INTO HMM_Hits(Hit_HASH, Protein_Accession, HMM_Model, HMM_Score, HMM_E_Value, Ali_From, Ali_To, HMM_From, HMM_To, HMM_Coverage) 
			              VALUES(?,?,?,?,?,?,?,?,?,?)''', hit)
#===========================================================================================================
# Main program code:
	
# House keeping...
argsCheck(5) # Checks if the number of arguments are correct.

# Stores file one for input checking.
print ">> Starting up..."
OrganismFile = sys.argv[1]
PhyloCSVFile = sys.argv[2]
HMMFile      = sys.argv[3]
sqlFile      = sys.argv[4]

HMMName = path.split(HMMFile)[1].rstrip(".hmm")
OrganismName = path.split(OrganismFile)[1].rstrip(".faa")

# File extension checks
print ">> Performing file extention checks..."
if not OrganismFile.endswith(".faa"):
	print "[Warning] " + OrganismFile + " may not be a fasta file!"
if not PhyloCSVFile.endswith(".csv"):
	print "[Warning] " + PhyloCSVFile + " may not be a csv file!"
if not HMMFile.endswith(".hmm"):
	print "[Warning] " + HMMFile + " may not be a HMM file!"
if not sqlFile.endswith(".sqlite"):
	print "[Warning] " + sqlFile + " may not be a sqlite file!"

# Read in genbank file as a sequence record object.
try:
	print ">> Opening FASTA File: " + OrganismFile
	handle = open(OrganismFile, "rU")
	records = SeqIO.parse(handle, "fasta")
	handle.close()
except IOError:
	print "Failed to open " + OrganismFile
	sys.exit(1)
	
# Opens OrganismDB CSV file for reading and stores as hash table.
OrganismHash = {}
try:
	print ">> Opening CSV File: " + PhyloCSVFile
	readFile = open(PhyloCSVFile, "r")
	reader  = csv.reader(readFile) # Opens file with csv module which takes into account verying csv formats and parses correctly.
	print ">> Good CSV file."
	for row in reader:
		OrganismHash[row[0]] = row
	readFile.close()
except IOError:
	print "Failed to open " + PhyloCSVFile
	sys.exit(1)

# Gets HMM Length.
try:
	HMMLengthRegex = re.compile("^LENG\s*\d*$")
	print ">> Opening HMM File: " + HMMFile
	with open(HMMFile, "rU") as inFile:
		currentLine = inFile.readline()
		while not HMMLengthRegex.match(currentLine):
			currentLine = inFile.readline()
		HMMLength = int(currentLine.split()[1])
except IOError:
	print "Failed to open " + HMMFile
	sys.exit(1)
	
print ">> Extracting Protein Annotations..."
AnnotationFASTADict = createProteomeHash(OrganismFile) # Creates a dictionary containing all protein annotations in the gbk file.

print ">> Extracting Organism Info..."
OrganismInfo = OrganismHash[OrganismName]

print ">> Running Hmmsearch..."
FASTAString = "".join(AnnotationFASTADict.values()) # Saves these annotations to a string.
HMMResults  = runHMMSearch(FASTAString, HMMFile) # Runs hmmsearch.

if (len(HMMResults) < 1):
	print "Hmmearch found no hits for " + HMMFile + " in " + OrganismName
	print "Aborting..."
	sys.exit(1)

print ">> Parsing and filtering hmmsearch results..."
HMMHitTable = parseHmmsearchResults(HMMResults, HMMName, HMMLength) # Parses hmmsearch results into a two dimensional array.
HMMHitTable = filterHMMHitTable(HMMHitTable)

if (len(HMMHitTable) < 1):
	print "We found no quality hits for " + HMMFile + " in " + OrganismName
	print "Aborting..."
	sys.exit(1)

print ">> Extracting Hit Proteins."
HitProteins = getHitProteins(HMMHitTable, AnnotationFASTADict, OrganismName) # Gets hit protein FASTAs.

if path.isfile(sqlFile):
	try:
		HMMDB = sqlite3.connect(sqlFile, timeout = 3600)
		print ">> Opened database successfully."
		cursor = HMMDB.cursor()
		
		print ">> Inserting Organism Info."
		insertOrganismInfo(cursor, OrganismInfo)
		print ">> Inserting Proteins."
		insertProteins(cursor, HitProteins)
		print ">> Inserting Hits."
		insertHits(cursor, HMMHitTable)
			
		HMMDB.commit()
		HMMDB.close()
	except sqlite3.Error as Error:
		print "sqlite3 Error: " + str(Error)
		print "The program will be aborted."
		sys.exit(1)		
else:
	print "Failed to open " + sqlFile
	sys.exit(1)
print ">> Done!"