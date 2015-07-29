#!/usr/bin/env python

"""
Created by: Lee Bergstrand

Description:    Functions for HMM search.
"""

# Imports & Setup:
import csv
import sys
from Bio import SeqIO
import subprocess
from multiprocessing import cpu_count
import re

processors = cpu_count()  # Gets number of processor cores for HMMER.

# Regex's
LocusRegex = re.compile("\(Locus:\s\S*\)")
LocationRegex = re.compile("\(Location:\s\[(\S*)\:(\S*)\]\((\S)\)\)")


# ----------------------------------------------------------------------------------------
def extract_sequence_records(organism_file_path, file_type):
	"""
	Read in sequence files as a sequence record object using Biopython.

	:param organism_file_path: The path to the input file.
	:return: Biopython sequence record object.
	"""
	try:
		print(">> Opening " + organism_file_path)
		handle = open(organism_file_path, "rU")
		try:
			records = list(SeqIO.parse(handle, file_type))
		except ValueError as error:
			print("Error has occurred while parsing " + organism_file_path + "!")
			print(str(error))
			sys.exit(1)
		handle.close()
	except IOError:
		print("Failed to open " + organism_file_path)
		sys.exit(1)
	return records


# ----------------------------------------------------------------------------------------
def generate_fasta_string(sec_record_list):
	"""
	Creates a FASTA formatted string containing sequences from a list of sequence record objects.

	:param sec_record_list: List of Biopython sequence record objects.
	:return: String containing FASTA formatted strings.
	"""
	fasta_string_list = []

	for record in sec_record_list:
		fasta_string_list.append(record.format("fasta"))

	fasta_string = ''.join(fasta_string_list)
	return fasta_string


# ----------------------------------------------------------------------------------------
def generate_fasta_dict(sec_record_list):
	"""
	Creates a dictionary containing FASTA formatted strings from a list of sequence record objects.
	This dictionary is keyed by the sequence ID.

	:param sec_record_list: List of Biopython sequence record objects.
	:return: Dictionary containing FASTA formatted strings.
	"""
	fasta_string_dict = {}

	for record in sec_record_list:
		fasta_string_dict[record.id] = record.format("fasta")

	return fasta_string_dict


# ----------------------------------------------------------------------------------------
def run_hmm_search(fasta_string, hmmer_model_path):
	"""
	Runs HMMER with settings specific for extracting subject sequences.

	:param fasta_string: String containing protein sequences in FASTA format.
	:param hmmer_model_path: Path to the HMM model to be used as a query.
	:return:
	"""
	process = subprocess.Popen(["hmmsearch", "--acc", "--cpu", str(processors), hmmer_model_path, "-"],
	                           stdin=subprocess.PIPE, stdout=subprocess.PIPE, bufsize=1)

	# This returns a list with both stderr and stdout. Only return stdout. Fail if error.
	stdout, error = process.communicate(fasta_string)
	if error:
		print(str(error))
		sys.exit(1)
	else:
		return stdout


# ------------------------------------------------------------------------------------------------------------
def get_hit_protein_data(hmm_hit_table, annotation_fasta_dict, organism_accession):
	"""
	Creates a list of lists which contain protein information.

	:param hmm_hit_table: Table of HMM hit objects.
	:param annotation_fasta_dict: Dictionary containing FASTA sequences keyed by their IDs
	:param organism_accession: Name of the organism
	:return:
	"""

	hit_proteins = []
	for hit in hmm_hit_table:
		protein_accession = hit.target_protein
		protein_fasta = annotation_fasta_dict[protein_accession]

		locus = str(LocusRegex.search(protein_fasta).group(0))
		locus = locus.split()[1].rstrip(")")
		location_data = LocationRegex.search(protein_fasta)
		try:
			start = int(location_data.group(1))
			end = int(location_data.group(2))
			strand = location_data.group(3)
			protein_data = [protein_accession, organism_accession, locus, start, end, strand, protein_fasta]
			hit_proteins.append(protein_data)
		except AttributeError as error:
			print(hit)
			print(protein_fasta)
			print(location_data)
			print("This is the organism: ", organism_accession)
			print("The AttributeError was ", str(error))
			sys.exit(1)
	return hit_proteins


def extract_csv_dict(input_csv_path):
	"""
	Opens OrganismDB CSV file for reading and stores as dictionary.

	:param input_csv_path: Path to the input OrganismDB CSV file.
	:return: Dictionary with each row in the CSV keyed by the organism accession (CSV row one).
	"""

	organism_hash = {}
	try:
		print(">> Opening CSV File: " + input_csv_path)
		read_file = open(input_csv_path, "r")
		reader = csv.reader(read_file)
		print(">> Good CSV file.")
		for row in reader:
			organism_hash[row[0]] = row  # Row[0] is the organism accession.
		read_file.close()
	except IOError:
		print("Failed to open " + input_csv_path)
		sys.exit(1)

	return organism_hash
