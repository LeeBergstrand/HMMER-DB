#!/bin/bash

accession=$(ls *.gbk | head -1 | sed "s/\.gbk//");

python /Users/lee/Dropbox/RandD/Repositories/HMMER-DB/GetOrganismTableRow.py $accession.gbk;
pigz -dc *.gbk.tgz | tar xf -;
find . -name "*.gbk" | parallel --gnu -j +0 python /Users/lee/Dropbox/RandD/Repositories/HMMER-DB/GetAnnotatedFASTA.py {};
cat *.faa > "$accession.faa"
mv ./OrganismDB.csv ./$accession.csv
mv $accession.faa $accession.csv ../