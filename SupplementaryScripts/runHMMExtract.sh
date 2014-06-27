#!/bin/bash
# A simple script for the batch running of HMMExtract.py.

timestamp=$(date +"%T")

for HMM;
do
	echo Running HMMExtract on $HMM
	time find ./ -name "*.faa" | parallel --gnu --progress -k -j +0 python ./SearchExtractLoad.py {} OrganismDB.csv $HMM ./TestDBs/TestDB.sqlite 2>&1 | tee -a "Log$timestamp.txt"
	date +"%T" >> "Log$timestamp.txt"
	echo 
done
echo All hmms searched.
 
exit 0

