#!/usr/bin/env python

"""
Created by: Lee Bergstrand

Description:    Parses HMM search output.

Requirements:   - This script requires HMMER 3.1 or later.
"""

# Imports & Setup:
import re
import sys
from hashlib import md5

hit_row_regex = re.compile("^\s*\d\s*((\?)|(\!))\s*")


# ========
# Classes:
# ========

class HMMHit(object):
	def __init__(self, target_protein, hmm_name, score, e_value, hmm_from, hmm_to, ali_from, ali_to, hmm_length):
		self.target_protein = str(target_protein)
		self.hmm_name = str(hmm_name)
		self.score = float(score)
		self.e_value = float(e_value)
		self.hmm_from = int(hmm_from)
		self.hmm_to = int(hmm_to)
		self.hmm_overlap = self.hmm_to - self.hmm_from
		self.ali_from = int(ali_from)
		self.ali_to = int(ali_to)
		self.ali_length = self.ali_to - self.ali_from
		self.hmm_coverage = float(self.hmm_overlap) / float(hmm_length)

	def __repr__(self):
		out_string = """
		===========================
		Target Protein: %s
		HMM Name:       %s
		Score:          %f
		e_value:        %f
		hmm_from:       %d
		hmm_to:         %d
		hmm_overlap     %d
		ali_from:       %d
		ali_to:         %d
		ali_length:     %d
		hmm_coverage:   %f
		""" % (
			self.target_protein, self.hmm_name, self.score, self.e_value, self.hmm_from, self.hmm_to,
			self.hmm_overlap, self.ali_from, self.ali_to, self.ali_length, self.hmm_coverage)

		return out_string

	def __str__(self):
		return self.__repr__()

	def get_md5(self):
		hash_string = "".join([str(x) for x in self.__dict__.values()])  # Join all attributes into a single string.
		hash_md5 = md5(hash_string.encode('utf-8')).hexdigest()  # Create md5 hash.
		return hash_md5


# ==========
# Functions:
# ==========

# ------------------------------------------------------------------------------------------------------------
def parse_hmmsearch_results(hmm_results_string, hmm_name, hmm_length):
	"""
	Parses HMM searches text output and generates a two dimensional array of the domain alignments results.

	:param hmm_results_string: hmmsearch text output as a string.
	:param hmm_name: The name of the HMM file.
	:param hmm_length: The length of the HMM file.
	:return: List of HMM hit objects.
	"""

	hmm_hit_list = hmm_results_string.split(">>")  # Splits output at domain alignments.
	del hmm_hit_list[0]  # Removed string content above first >> which does not need to be iterated through.

	# Removes alignment visualization by spiting on the visualization header
	# and keeping all that is above it (first element) which is the HMM domain table.
	hmm_hit_list_cleaned = [x.split("Alignments")[0] for x in hmm_hit_list]

	hmm_object_list = []
	for protein in hmm_hit_list_cleaned:
		domain_table = protein.splitlines()
		target_protein_name = domain_table[0].split()[0]  # Gets target protein accession from domain table header row.

		for row in domain_table:
			if hit_row_regex.match(row):  # If row is a domain table line.
				column = row.split()

				hmm_hit_object = HMMHit(
					target_protein=target_protein_name,
					hmm_name=hmm_name,
					score=column[2],
					e_value=column[5],
					hmm_from=column[6],
					hmm_to=column[7],
					ali_from=column[9],
					ali_to=column[10],
					hmm_length=hmm_length)

				hmm_object_list.append(hmm_hit_object)

	return hmm_object_list


# ------------------------------------------------------------------------------------------------------------
def filter_hmm_hit_list(hmm_hit_list, e_value_cutoff="1e-25", hmm_coverage=0.3, max_align_overlap=0.5):
	"""
	Filters HMM gits by E-Value, Coverage and Overlap between hits.

	:param hmm_hit_list: List of HMM hit objects.
	:param e_value_cutoff: The E-Value cutoff for hits.
	:param hmm_coverage: The HMM coverage cutoff for hits.
	:param max_align_overlap: The maximum overlap percentage between overlapping HMM hits.
	:return: List of filtered HMM hit objects.
	"""
	hmm_hit_list = [hit for hit in hmm_hit_list if hit.e_value < float(e_value_cutoff)]  # Filters hits by E-value.

	i = 0
	while i < (len(hmm_hit_list) - 1):
		hit_one = hmm_hit_list[i]  # Current Row in hit table.
		hit_two = hmm_hit_list[i + 1]  # Row below.
		if hit_one.target_protein == hit_two.target_protein:
			overlap_between_hits = hit_one.ali_to - hit_two.ali_from
			if overlap_between_hits > 0:
				# If the overlap is greater than 50% of either alignment.
				if ((float(overlap_between_hits) / float(hit_one.ali_length)) > max_align_overlap) or (
							(float(overlap_between_hits) / float(hit_two.ali_length)) > max_align_overlap):

					if hit_one.e_value < hit_two.e_value:
						hmm_hit_list.remove(hit_two)
					else:
						hmm_hit_list.remove(hit_one)
					i -= 1  # Resets list index.
		i += 1

	hmm_hit_list = [hit for hit in hmm_hit_list if hit.hmm_coverage > hmm_coverage]  # Filters by Query Coverage.
	return hmm_hit_list


# ------------------------------------------------------------------------------------------------------------
def get_hmm_length(hmm_path):
	"""
	Gets HMM Length.

	:param hmm_path: path to the HMM file.
	:return: The length of the hmm.
	"""
	try:
		hmm_length_regex = re.compile("^LENG\s*\d*$")
		print(">> Opening HMM file: " + hmm_path)
		with open(hmm_path, "rU") as inFile:
			current_line = inFile.readline()
			while not hmm_length_regex.match(current_line):
				current_line = inFile.readline()
			hmm_length = int(current_line.split()[1])
	except IOError:
		print("Failed to open " + hmm_path)
		sys.exit(1)

	return hmm_length
