#!/usr/bin/env python

"""
Created by: Lee Bergstrand

Description:    A program that extracts the protein annotations from a fasta file and searches these
				annotations using HMMsearch and an HMM file. It then stores hits along with organism
				information (gathered from a csv file) in a sqlite3 database.

Requirements:   - This script requires the Biopython module: http://biopython.org/wiki/Download
				- This script requires HMMER 3.1 or later.
				- This script requires sqlLite3 or later. HMMExtract.py <organism.faa> <organisms.csv> <hmm.hmm> <sqldb.sqlite>
"""

# Imports & Setup:
import argparse
import sqlite3
from hmm_parser import *
from lib import *
from os import path
from multiprocessing import cpu_count

processors = cpu_count()  # Gets number of processor cores for HMMER.


# ==========
# Functions:
# ==========


# ----------------------------------------------------------------------------------------
def main(args):
	input_file_path = args.in_file[0]
	organism_csv_path = args.in_csv[0]
	hmm_model_paths = args.hmm_models
	database_path = args.database[0]

	print('\nHMMER-DB')
	print('====================================================')

	check_extensions(input_file_path, organism_csv_path, hmm_model_paths, database_path)

	print('')

	sequence_record_list = extract_sequence_records(input_file_path, 'fasta')
	fasta_string = generate_fasta_string(sequence_record_list)
	fasta_dict = generate_fasta_dict(sequence_record_list)
	organism_data_dict = extract_csv_dict(organism_csv_path)
	organism_file_name = path.basename(input_file_path).split('.')[0]
	organism_data = organism_data_dict[organism_file_name]

	organism_accession = organism_data[0]

	print('')

	hmm_hit_list = []
	protein_data_list = []
	for hmm_path in hmm_model_paths:
		hits_to_add, proteins_to_add = run_hmm_search(hmm_path, fasta_string, fasta_dict,
		                                              organism_accession, processors)

		hmm_hit_list.extend(hits_to_add)
		protein_data_list.extend(proteins_to_add)

	# To account for sqlite3 table lock on write timeout
	# is delayed in proportion to the number of CPUs used.
	timeout_for_parallelism = 225 * processors
	if path.isfile(database_path):
		try:
			print('')

			hmm_db = sqlite3.connect(database_path, timeout=timeout_for_parallelism)
			print(">> Opening sqlite3 file: " + database_path)
			cursor = hmm_db.cursor()

			print(">> Inserting organism info...")
			insert_organism_info(cursor, organism_data)
			print(">> Inserting cached protein data...")
			insert_proteins(cursor, protein_data_list)
			print(">> Inserting cached HMM hit data...")
			insert_hits(cursor, hmm_hit_list)

			hmm_db.commit()
			hmm_db.close()
		except sqlite3.Error as error:
			print("sqlite3 Error: " + str(error))
			print("The program will be aborted.")
			sys.exit(1)
	else:
		print("Failed to open " + database_path)
		sys.exit(1)
	print("\n>> Done!")


# ------------------------------------------------------------------------------------
def run_hmm_search(hmm_path, fasta_string, fasta_dict, organism_accession, processes):
	"""
	Runs the HMM search using hmmsearch.

	:param hmm_path: Path to the HMM search file.
	:param fasta_string: A FASTA string containing proteins to be searched.
	:param fasta_dict: A dict of FASTA strings containing proteins to be searched keyed by sequence ID.
	:param organism_accession: The accession of the organism.
	:param processes: The number of threads to use when running hmmsearch.
	:return: list of HMM hit objects and list of lists of protein data.
	"""

	print('')

	hmm_name = path.basename(hmm_path).split('.')[0]
	hmm_length = get_hmm_length(hmm_path)

	print('>> Running hmmsearch on ' + str(processes) + ' CPUs...')
	hmm_results_string = hmm_search(fasta_string, hmm_path, processes)
	hmm_hit_list = parse_hmmsearch_results(hmm_results_string, hmm_name, hmm_length)

	print('>> Filtering HMM hits...')
	filtered_hmm_hit_list = filter_hmm_hit_list(hmm_hit_list)

	protein_data_list = get_hit_protein_data(filtered_hmm_hit_list, fasta_dict, organism_accession)

	print('>> Caching HMM hit and subject proteins data...')

	return filtered_hmm_hit_list, protein_data_list


# ----------------------------------------------------------------------------------------
if __name__ == '__main__':
	descriptor = """
	A program that extracts the protein annotations from a fasta file and searches these
	annotations using HMMsearch and an HMM file. It then stores hits along with organism
	information (gathered from a csv file) in a sqlite3 database.
	"""
	parser = argparse.ArgumentParser(description=descriptor)

	parser.add_argument('-i', '--in_file', metavar='FASTA', nargs=1, help='''
	The input FASTA file containing protein sequences (Created by GenbankToFASTAandOrganismTableRow.py).''')

	parser.add_argument('-c', '--in_csv', metavar='CSV', nargs=1, help='''
	The CSV file containing the information of all input organism (Created by GenbankToFASTAandOrganismTableRow.py).''')

	parser.add_argument('-m', '--hmm_models', metavar='HMM', nargs='+', help='''
	The HMM model files representing proteins''')

	parser.add_argument('-d', '--database', metavar='DATABASE', nargs=1, help='''
	The input sqlite3 database for which the organism info and HMM results are writen to.''')

	parser.add_argument('-p', '--processes', metavar='PROCESSES', nargs=1, default=processors, help='''
	Number of parallel processes to be used by hmmsearch.''')

	cli_args = parser.parse_args()

	# At minimum we require all CLI inputs.
	proceed = True

	if cli_args.in_file is None:
		print("Error: Missing input FASTA file path...")
		proceed = False

	if cli_args.in_csv is None:
		print("Error: Missing CSV file path...")
		proceed = False

	if cli_args.hmm_models is None:
		print("Error: Missing HMM file paths...")
		proceed = False

	if cli_args.database is None:
		print("Error: Missing sqlite3 database path...")
		proceed = False

	if proceed:
		main(cli_args)
	else:
		print("")
		parser.print_help()
		print("")
