#!/usr/bin/env python 
# Created by: Lee Bergstrand 
# Descript: A script that extracts the FASTA sequences of proteins with hits for specific Hidden Markov Model
#           within a HMMER-DB database. 
#
# Requirements: - This script requires the Biopython module: http://biopython.org/wiki/Download
#               - This script requires sqlLite3 or later.
#
# Usage: FastaExtract.py <HMMName> <sqldb.sqlite>
# Example: FastaExtract.py <HsaCGramPos> steriodDB.sqlite
#----------------------------------------------------------------------------------------

import sys
import sqlite3

