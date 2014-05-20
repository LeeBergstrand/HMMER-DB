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
import tempfile
import cStringIO
from os import path
from Bio import SeqIO
from Bio import SearchIO  
from multiprocessing import cpu_count

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
def runHMMSearch(FASTA, HMMERDBFile, tempHMMTTabFile):
	Found16S = True
	process = subprocess.Popen(["hmmsearch", "--acc", "--cpu", str(processors), "--pfamtblout", tempHMMTTabFile.name, HMMERDBFile, "-"], stdin = subprocess.PIPE, stdout = subprocess.PIPE, bufsize = 1)
	stdout, error = process.communicate(FASTA) # This returns a list with both stderr and stdout. Only want stdout which is first element.
	if error:
		print error
		sys.exit(1)
	#else:
		#print stdout
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
#===========================================================================================================
# Main program code:
	
# House keeping...
argsCheck(4) # Checks if the number of arguments are correct.

# Stores file one for input checking.
print ">> Starting up..."
OrganismFile = sys.argv[1]
HMMFile      = sys.argv[2]
sqlFile      = sys.argv[3]

dirname, basename = path.split(OrganismFile)
basename = basename.rstrip(".gbk")

tempHMMTTabFile = tempfile.NamedTemporaryFile(prefix = basename, dir = dirname, suffix = ".hmmer") # Creates temporary for hmmer to write to.

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

print ">> Extracting Protein Annotations..."
AnnotationFASTADict = getProtienAnnotationFasta(record) # Creates a dictionary containing all protein annotations in the gbk file.
FASTAString = "".join(AnnotationFASTADict.values()) # Saves these annotations to a string.
open("Blam.faa", "w").write(FASTAString)
runHMMSearch(FASTAString, HMMFile, tempHMMTTabFile) # Runs hmmer and writes to temporary file.

HMMResults = tempHMMTTabFile.read()
#print HMMResults
HMMResults = HMMResults.split("Domain scores")[1] # Extracts domain scores.
#print HMMResults
HMMResults = HMMResults.splitlines()

HMMResultsCleaned = []
for i in HMMResults:
	if "#" not in i:
		HMMResultsCleaned.append(i)
del HMMResultsCleaned[0] # Deletes empty element at start of list.
	
HMMResultsCleaned = [x.split() for x in HMMResultsCleaned]
	
# print HMMResultsCleaned

