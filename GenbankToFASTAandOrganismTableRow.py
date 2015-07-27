#!/usr/bin/env python

"""
Created by: Lee Bergstrand

Description:    A program that extracts the proteins annotations from a Genbank file and as well as some
				information about the organism in the file. Stores the protein annotations as a Fasta.
				Appends a csv file with the organism info.

Requirements: - This script requires the Biopython module: http://biopython.org/wiki/Download
"""

# Imports & Setup:
import sys
from Bio import SeqIO
import os
import argparse


# ==========
# Functions:
# ==========

# ----------------------------------------------------------------------------------------
def main(args):
	input_file_path = args.in_file[0]
	organism_info_flag = args.organism_info
	annotation_flag = args.annotation

	fasta_file_name = os.path.join(os.getcwd(), str(os.path.basename(input_file_path).split('.')[0]) + '.faa')
	csv_file_name = os.path.join(os.getcwd(), 'OrganismDB.csv')

	check_extension(input_file_path)
	seq_record = extract_sequence_record(input_file_path)

	if organism_info_flag:
		csv_row = get_organism_info(seq_record)
		write_csv(csv_row, csv_file_name)

	if annotation_flag:
		fasta = get_coding_annotation_fasta(seq_record)
		write_fasta(fasta, fasta_file_name)


# ----------------------------------------------------------------------------------------
def check_extension(in_path):
	"""
	Checks input file extension and give a warning if it is incorrect.

	:param in_path: The path to the input file.
	"""
	in_file_extension = os.path.splitext(in_path)[-1]
	if not in_file_extension == ".gbk":
		print("[Warning] " + in_file_extension + " may not be a Genbank file!")


# ----------------------------------------------------------------------------------------
def extract_sequence_record(organism_file):
	"""
	Read in Genbank file as a sequence record object using Biopython.

	:param organism_file: The path to the input file.
	:return: Biopython sequence record object.
	"""
	try:
		print(">> Opening " + organism_file)
		handle = open(organism_file, "rU")
		try:
			record = SeqIO.read(handle, "genbank")
		except ValueError as error:
			print("Error has occurred while parsing " + organism_file + "!")
			print(str(error))
			sys.exit(1)
		handle.close()
	except IOError:
		print("Failed to open " + organism_file)
		sys.exit(1)
	return record


# ------------------------------------------------------------------------------------------------------------
def get_organism_info(seq_record):
	"""
	When passed a sequence record object returns a string containing the organisms info.

	:param seq_record: A Biopython sequence record object.
	:return: A string containing a CSV file row.
	"""
	print(">> Extracting Organism Info...")
	organism_id = seq_record.id
	description = seq_record.description.replace(",", "")
	if 'plasmid' in description.lower():
		accession_type = 'Plasmid'
	elif 'chromosome' in description.lower():
		accession_type = 'Chromosome'
	elif 'genome' in description.lower():
		accession_type = 'Chromosome'
	else:
		accession_type = 'Unknown'
	source = seq_record.annotations['source'].replace(",", "")
	taxonomy = "_".join(seq_record.annotations['taxonomy'])
	organism_genome_length = len(seq_record.seq)

	organism_string = organism_id + "," + accession_type + "," + description + "," + source + "," + taxonomy + "," + str(
		organism_genome_length) + "\n"
	return organism_string


# -------------------------------------------------------------------------------------------------------------
def get_coding_annotation_fasta(seq_record):
	"""
	When passed a sequence record object returns an array of FASTA strings for each annotation.

	:param seq_record: A Biopython sequence record object.
	:return: A FASTA file string containing all sequences record object CDS sequence features.
	"""
	fasta = []
	features = seq_record.features  # Each sequence has a list (called features) that stores seqFeature objects.
	for feature in features:  # For each feature on the sequence
		if feature.type == "CDS":  # CDS means coding sequence (These are the only feature we're interested in)
			feat_qualifiers = feature.qualifiers  # Each feature contains a dictionary called qualifiers which contains
			# data about the sequence feature (for example the translation)

			start = int(feature.location.start)  # Type-casting to int strips fuzzy < > characters.
			end = int(feature.location.end)
			strand = feature.location.strand
			if strand is None:
				strand = "?"
			elif int(strand) < 0:
				strand = "-"
			elif int(strand) > 0:
				strand = "+"
			else:
				strand = "?"

			location = "[" + str(start) + ":" + str(end) + "](" + strand + ")"

			# Gets the required qualifiers. Uses featQualifiers.get to return the qualifiers or a default value if the qualifiers
			# is not found. Calls strip to remove unwanted brackets and ' from qualifiers before storing it as a string.
			protein_id = str(feat_qualifiers.get('protein_id', 'no_protein_id')).strip('\'[]')
			if protein_id == 'no_protein_id':
				continue  # Skips the iteration if protein has no id.
			protein_locus = str(feat_qualifiers.get('locus_tag', 'no_locus_tag')).strip('\'[]')
			gene = str(feat_qualifiers.get('gene', 'no_gene_name')).strip('\'[]')
			product = str(feat_qualifiers.get('product', 'no_product_name')).strip('\'[]')
			translated_protein = str(feat_qualifiers.get('translation', 'no_translation')).strip('\'[]')

			fasta_part_one = ">" + protein_id + " " + gene + "-" + product + " (Locus: " + protein_locus + ")"
			fasta_part_two = " (Location: " + location + ")" + "\n" + translated_protein + "\n"
			fasta.append(fasta_part_one + fasta_part_two)
	fasta_string = "".join(fasta)
	return fasta_string


# ----------------------------------------------------------------------------------------
def write_fasta(fasta, out_file_path):
	"""
	Writes annotations to FASTA file.

	:param fasta:
	:param out_file_path:
	"""
	try:
		print(">> Writing fasta File...")
		fasta_writer = open(out_file_path, "w")
		fasta_writer.write(fasta)
		fasta_writer.close()
	except IOError:
		print("Failed to open " + out_file_path)
		sys.exit(1)


# ----------------------------------------------------------------------------------------
def write_csv(organism_string, out_file_path):
	"""
	Appends to organism CSV.

	:param organism_string:
	:param out_file_path:
	"""
	try:
		print(">> Writing to organism info to CSV file...")
		writer = open(out_file_path, "a")
		writer.write(organism_string)
		writer.close()
	except IOError:
		print("Failed to open " + "OrganismDB.csv")
		sys.exit(1)
	print(">> Done")


# ----------------------------------------------------------------------------------------
if __name__ == '__main__':
	descriptor = """
	A program that extracts the proteins annotations from a Genbank file and as well as some
	information about the organism in the file. Stores the protein annotations as a Fasta.
	Appends a csv file with the organism info.
	"""

	parser = argparse.ArgumentParser(description=descriptor)

	parser.add_argument('-i', '--in_file', metavar='GENBANK', nargs=1, help='''
	The input Genbank file.''')

	parser.add_argument('-O', '--organism_info', action='store_true', help='''
	Flag to write an organism table row.''')

	parser.add_argument('-A', '--annotation', action='store_true', help='''
	Flag  to extract CDS from input file and write them to a FASTA file.''')

	cli_args = parser.parse_args()

	# At minimum we require an input file, and one CLI flag.
	proceed = True

	if cli_args.in_file is None:
		print("Error: Missing sequence list path...")
		proceed = False

	if not cli_args.organism_info and not cli_args.annotation:
		print("Error: Missing sequence FASTA file path...")
		print("Error: Missing query BLAST database path...")
		proceed = False

	if proceed:
		main(cli_args)
	else:
		print("")
		parser.print_help()
		print("")
