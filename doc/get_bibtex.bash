#!/usr/bin/env bash

########################################################
#
#    Copyright (c) 2015,2019-2020,2022,2025
#      SMASH Team
#
#    GNU General Public License (GPLv3 or later)
#
########################################################

# This script produces a bibtex file 'inspire.bib' including all papers
# referenced in the Smash source code via the doxygen command \iref{}.
#
# It first extracts all bibtex keys from the source and then queries Inspire
# for each of them, in order to get the full bibtex entry.

# Set up a trap on exit to clean up temporary files
trap 'rm -f ${tmp_bib_file}; echo' EXIT

# Make the script callable from anywhere (by changing to script directory)
cd "$(dirname "${BASH_SOURCE[0]}")"; echo ''

# A couple global variable
readonly tmp_bib_file='tmp.bib'
readonly bib_file='inspire.bib'

# Extract keys cited in documentation
readonly cited_keys=( $(grep -hor '\\iref{[A-Za-z0-9:-]*}' ../src | sort -u | sed 's/\\iref{\([^}]*\)}/\1/g') )
readonly number_of_keys=${#cited_keys[@]}

# Iterate over all keys and extract all bibtex entries
fetched_keys=0
failed_to_fetch=0
unfetched_keys=()
for key in "${cited_keys[@]}"; do
  printf "\n\e[1A\e[K Fetching key %${#number_of_keys}d/%d [%d error(s) occurrred] -> key = %s"\
         $((fetched_keys+failed_to_fetch+1)) ${number_of_keys} ${failed_to_fetch} "${key}"
  # Search for bibtex key on Inspire, obtain it using the API
  #  -> https://github.com/inspirehep/rest-api-doc#internal-identifiers
  bibtex_entry=$(wget -q -O - "https://inspirehep.net/api/literature?q=${key}&format=bibtex")
  if [[ "${bibtex_entry}" != '' ]]; then
    echo "${bibtex_entry}" >> "${tmp_bib_file}"
    ((fetched_keys++))
  else
    unfetched_keys+=( "${key}")
    ((failed_to_fetch++))
  fi
done

# Print report and move temporary in place if successful
printf '\n\e[1A\e[K Successfully fetched %d/%d keys\n' ${fetched_keys} ${#cited_keys[@]}
if [[ ${failed_to_fetch} -eq 0 ]]; then
  mv "${tmp_bib_file}" "${bib_file}"
else
  printf '\n Unable to fetch following keys:\n'
  printf '  - %s\n' "${unfetched_keys[@]}"
fi
