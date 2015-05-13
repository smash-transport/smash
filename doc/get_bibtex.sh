#!/bin/bash
#
# This script produces a bibtex file 'smash.bib' including all papers
# referenced in the Smash source code via the doxygen command \iref{}.
#
# It first extracts all bibtex keys from the source and then queries Inspire
# for each of them, in order to get the full bibtex entry.

# make sure no earlier bib file is present
rm -f smash.bib

# extract all bibtex keys from the source files
for line in `grep -hor '\iref{.*}' ../src | sort -u | sed 's/iref{//' | sed 's/}//'`; do
  echo "fetching $line ..."
  # search for bibtex key on Inspire
  wget -q -O tmp.txt "http://inspirehep.net/search?p=$line&of=hx"
  # extract bib entry add append it to bibtex file
  sed -n "/<pre>/,/<\/pre>/p" tmp.txt | head -n -1 | tail -n +2 >> smash.bib
done

# clean up
rm tmp.txt
