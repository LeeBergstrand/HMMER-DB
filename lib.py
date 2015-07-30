#!/usr/bin/env python

"""
Created by: Lee Bergstrand

Description:    Functions for HMMER-DB.
"""

# Imports & Setup:
import csv
import sys
from Bio import SeqIO
import subprocess
import re

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
		print(">> Opening FASTA file: " + organism_file_path)
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


# -----------------------------------------------------------------------------------------------------------
def check_extensions(organism_file_path, csv_file_path, hmm_file_paths, sql_file_paths):
	"""
	Performs file extension checks.

	:param organism_file_path: Path to the organism database file.
	:param csv_file_path: Path to the organism information database file.
	:param hmm_file_paths: Path to the HMM model file.
	:param sql_file_paths: Path to the sqlite3 file.
	"""

	print(">> Performing file extension checks...")
	if not organism_file_path.endswith(".faa"):
		print("[Warning] " + organism_file_path + " may not be a fasta file!")
	if not csv_file_path.endswith(".csv"):
		print("[Warning] " + csv_file_path + " may not be a csv file!")

	for hmm_path in hmm_file_paths:
		if not hmm_path.endswith(".hmm"):
			print("[Warning] " + hmm_path + " may not be a HMM file!")

	if not sql_file_paths.endswith(".sqlite"):
		print("[Warning] " + sql_file_paths + " may not be a sqlite file!")


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
def hmm_search(fasta_string, hmmer_model_path, processes):
	"""
	Runs HMMER with settings specific for extracting subject sequences.

	:param fasta_string: String containing protein sequences in FASTA format.
	:param hmmer_model_path: Path to the HMM model to be used as a query.
	:return: String containing hmmsearch output.
	"""
	process = subprocess.Popen(["hmmsearch", "--acc", "--cpu", str(processes), hmmer_model_path, "-"],
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
	:param organism_accession: The accession of the organism.
	:return: A list of lists of hit protein properties.
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


# -----------------------------------------------------------------------------------------------------------
def extract_csv_dict(input_csv_path):
	"""
	Opens OrganismDB CSV file for reading and stores as dictionary.

	:param input_csv_path: Path to the input OrganismDB CSV file.
	:return: Dictionary with each row in the CSV keyed by the organism accession (CSV row one).
	"""

	organism_data_csv = {}
	try:
		print(">> Opening organism CSV file: " + input_csv_path)
		read_file = open(input_csv_path, "r")
		reader = csv.reader(read_file)
		for row in reader:
			organism_data_csv[row[0].split('.')[0]] = row  # Row[0] is the organism accession.
		read_file.close()
	except IOError:
		print("Failed to open " + input_csv_path)
		sys.exit(1)

	return organism_data_csv


# -----------------------------------------------------------------------------------------------------------
def insert_organism_info(db_cursor, organism_info):
	"""
	Inserts organism info into DB.

	:param db_cursor: Sqlite3 database cursor.
	:param organism_info: List containing organism info.
	"""

	query = '''INSERT OR REPLACE INTO Organisms
(
	Organism_Accession,
	Accession_Type,
	Organism_Description,
	Source,
	Organism_Phylogeny,
	Sequence_Length
)
VALUES
	(?,?,?,?,?,?)'''

	db_cursor.execute(query, organism_info)


# -----------------------------------------------------------------------------------------------------------
# 8:
def insert_proteins(db_cursor, hit_proteins):
	"""
	Inserts protein info into DB.

	:param db_cursor: Sqlite3 database cursor.
	:param hit_proteins: List containing protein info.
	"""

	query = '''INSERT OR REPLACE INTO Proteins
(
	Protein_Accession,
	Organism_Accession,
	Locus,
	Start,
	"End",
	Strand,
	FASTA_Sequence
)
VALUES
	(?,?,?,?,?,?,?)'''

	for protein in hit_proteins:
		db_cursor.execute(query, protein)


# -----------------------------------------------------------------------------------------------------------
def insert_hits(cursor, hmm_hit_list):
	"""
	Inserts hits into DB and creates md5 hash for primary key.

	:param cursor: Sqlite3 database cursor.
	:param hmm_hit_list: List of hmm hit objects.
	"""

	query = '''INSERT OR REPLACE INTO HMM_Hits
(
	Hit_HASH,
	Protein_Accession,
	HMM_Model,
	HMM_Score,
	HMM_E_Value,
	Ali_From,
	Ali_To,
	HMM_From,
	HMM_To,
	HMM_Coverage
)
VALUES
	(?,?,?,?,?,?,?,?,?,?)'''

	for hit in hmm_hit_list:
		hit_list = [hit.get_md5(), hit.target_protein, hit.hmm_name, hit.score, hit.e_value, hit.ali_from, hit.ali_to,
		            hit.hmm_from, hit.hmm_to, hit.hmm_coverage]
		cursor.execute(query, hit_list)
