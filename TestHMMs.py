#!/usr/bin/env python

"""
Created by: Lee Bergstrand

Description:    Tests that HMMs finds all reference proteins in a given organism.

Requirements:   - This software requires the Biopython module: http://biopython.org/wiki/Download
				- This software requires HMMER 3.1 or later.
				- This software requires sqlLite3 or later.
"""

# Imports & Setup:
import argparse
from hmm_parser import *
from lib import *
from os import path
from multiprocessing import cpu_count

cpu_num = cpu_count()  # Gets number of processor cores for HMMER.


# ==========
# Functions:
# ==========


# ----------------------------------------------------------------------------------------
def main(args, processors):
	input_file_path = str(args.in_file[0])
	reference_accessions = str(args.list[0])
	hmm_model_paths = list(args.hmm_models)
	user_processes = int(args.processes[0])

	"""
	Only use the user specified process count if it is less
	than the number of cpu cores. Use all cpu cores by default.
	"""
	if 0 < user_processes < processors:
		processors = user_processes

	print("\nTesting HMMs")
	print('=============================================')

	sequence_record_list = extract_sequence_records(input_file_path, 'fasta')
	fasta_string = generate_fasta_string(sequence_record_list)

	hmm_hit_accessions = []
	for hmm_path in hmm_model_paths:
		print('')

		hmm_name = path.basename(hmm_path).split('.')[0]
		hmm_length = get_hmm_length(hmm_path)

		print('>> Running hmmsearch on ' + str(processors) + ' CPUs...')
		hmm_results_string = hmm_search(fasta_string, hmm_path, processors)
		hmm_hit_list = parse_hmmsearch_results(hmm_results_string, hmm_name, hmm_length)

		print('>> Filtering HMM hits...')
		filtered_hmm_hit_list = filter_hmm_hit_list(hmm_hit_list)

		print('>> Caching HMM hits...')
		hmm_hit_accessions.extend([hmm_hit.target_protein for hmm_hit in filtered_hmm_hit_list])

	hmm_hit_accessions = uniq_list(hmm_hit_accessions)
	reference_accessions = get_reference_accessions(reference_accessions)

	missing_accessions = []
	for reference in reference_accessions:
		if reference not in hmm_hit_accessions:
			missing_accessions.append(reference)

	print('')
	if not missing_accessions:
		print("No reference proteins are missing.")
	else:
		print("The following reference proteins are missing:")
		print('=============================================\n')
		for protein in missing_accessions:
			print('>' + protein)


# ----------------------------------------------------------------------------------------
def uniq_list(in_list):
	"""
	Takes a list of elements and removes duplicates and returns this list

	:param in_list: Input list
	:return: List containing unique items.
	"""

	seen = set()
	result = []
	for item in in_list:
		if item in seen:
			continue
		seen.add(item)
		result.append(item)
	return result


# ----------------------------------------------------------------------------------------
def get_reference_accessions(ref_file_path):
	"""
	Parses file containing a list of reference protein accessions.

	:param ref_file_path: Path to the reference protein file.
	:return: List of reference accessions.
	"""

	try:
		print(">> Opening organism CSV file: " + ref_file_path)
		read_file = open(ref_file_path, "r")
		data = read_file.read().splitlines()
		read_file.close()
	except IOError:
		print("Failed to open " + ref_file_path)
		sys.exit(1)
	return data


# ----------------------------------------------------------------------------------------
if __name__ == '__main__':
	descriptor = """
	Checking if a set of HMMs find all reference proteins found in an organism.
	"""
	parser = argparse.ArgumentParser(description=descriptor)

	parser.add_argument('-i', '--in_file', metavar='FASTA', nargs=1, help='''
	The input FASTA file containing protein sequences (Created by GenbankToFASTAandOrganismTableRow.py).''')

	parser.add_argument('-l', '--list', metavar='LIST', nargs=1, help='''
	List of required ''')

	parser.add_argument('-m', '--hmm_models', metavar='HMM', nargs='+', help='''
	The HMM model files representing proteins''')

	parser.add_argument('-p', '--processes', metavar='PROCESSES', nargs=1, default=[cpu_num], help='''
	Number of parallel processes to be used by hmmsearch.''')

	cli_args = parser.parse_args()

	# At minimum we require all CLI inputs.
	proceed = True

	if cli_args.in_file is None:
		print("Error: Missing input FASTA file path...")
		proceed = False

	if cli_args.list is None:
		print("Error: Missing accession list file path...")
		proceed = False

	if cli_args.hmm_models is None:
		print("Error: Missing HMM file paths...")
		proceed = False

	if proceed:
		main(cli_args, cpu_num)
	else:
		print("")
		parser.print_help()
		print("")
