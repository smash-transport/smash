#!/bin/bash
#
# This script produces a bibtex file 'smash.bib' including all papers
# referenced in the Smash source code via the doxygen command \iref{}.
#
# It first extracts all bibtex keys from the source and then queries Inspire
# for each of them, in order to get the full bibtex entry.

# make the script callable from anywhere (by changing to script directory)
cd `dirname \`which "$0"\``

# make sure no earlier bib file is present
rm -f smash.bib

# extract all bibtex keys from the source files
for line in `grep -hor '\iref{.*}' ../src | sort -u | sed 's/iref{//; s/}//'`; do
  echo "fetching $line ..."
  # search for bibtex key on Inspire, extract bib entry and append it to bibtex file
  wget -q -O - "http://inspirehep.net/search?p=$line&of=hx" | \
      sed -n "/<pre>/,/<\/pre>/p" | sed '$ d; 1d' >> smash.bib
done
