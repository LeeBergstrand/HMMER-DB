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
from hmm_parser import *
from lib import *
from os import path

# ==========
# Functions:
# ==========


# ----------------------------------------------------------------------------------------
def main(args):
	input_file_path = args.in_file[0]
	organism_csv_path = args.in_csv[0]
	hmm_model_path = args.hmm_model[0]
	database_path = args.database[0]

	sequence_record_list = extract_sequence_records(input_file_path, 'fasta')
	fasta_string = generate_fasta_string(sequence_record_list)
	fasta_dict = generate_fasta_dict(sequence_record_list)

	organism_accession = path.basename(input_file_path).split('.')[0]
	hmm_name = path.basename(hmm_model_path).split('.')[0]
	hmm_length = get_hmm_length(hmm_model_path)

	hmm_results_string = run_hmm_search(fasta_string, hmm_model_path)

	hmm_hit_list = parse_hmmsearch_results(hmm_results_string, hmm_name, hmm_length)
	filtered_hmm_hit_list = filter_hmm_hit_list(hmm_hit_list)

	protein_data = get_hit_protein_data(filtered_hmm_hit_list, fasta_dict, organism_accession)

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

	parser.add_argument('-m', '--hmm_model', metavar='HMM', nargs=1, help='''
	The HMM file representing a specific protein''')

	parser.add_argument('-d', '--database', metavar='DATABASE', nargs=1, help='''
	The input sqlite3 database for which the organism info and HMM results are writen to.''')

	cli_args = parser.parse_args()

	# At minimum we require all CLI inputs.
	proceed = True

	if cli_args.in_file is None:
		print("Error: Missing input FASTA file path...")
		proceed = False

	if cli_args.in_csv is None:
		print("Error: Missing CSV file path...")
		proceed = False

	if cli_args.hmm_model is None:
		print("Error: Missing HMM file path...")
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
