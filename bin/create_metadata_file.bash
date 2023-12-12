#!/usr/bin/env bash

#===================================================
#
#    Copyright (c) 2023
#      SMASH Team
#
#    GNU General Public License (GPLv3 or later)
#
#===================================================

trap 'printf "\n"' EXIT

#===================================================================
# This script expects a ".authors.json" file to be present were run
# and produces a "".zenodo.json" and/or a "AUTHORS.md" file in the
# same place. The jq tool is used to deal with JSON files and the
# yq tool is used to deal with YAML files.
#===================================================================

# Bash stricter mode
set -euo pipefail
shopt -s inherit_errexit

readonly \
    INPUT_FILE='.authors.json'\
    OUTPUT_AUTHORS_FILE='AUTHORS.md'\
    OUTPUT_ZENODO_FILE='.zenodo.json'\
    OUTPUT_CITATION_FILE='CITATION.cff'\
    SMASH_VERSION='3.0'\
    SMASH_RELEASE_DATE='2023-04-27'

CREATE_ZENODO_FILE='FALSE'
CREATE_AUTHORS_FILE='FALSE'
CREATE_CITATION_FILE='FALSE'

function main()
{
    parse_command_line_arguments "$@"
    make_preliminary_checks
    create_authors_file_if_requested
    create_zenodo_file_if_requested
    create_citation_file_if_requested
}

function make_preliminary_checks()
{
    if ! hash jq; then
        fail\
            "Program 'jq' not found."\
            'Follow instructions at https://jqlang.github.io/jq/download/ to install it.'
    fi
    if [[ ${CREATE_CITATION_FILE} = 'TRUE' ]]; then
        if ! hash yq; then
            fail\
                "Program 'yq' not found."\
                'Follow instructions at https://github.com/mikefarah/yq#install to install it.'
        fi
    fi
    if [[ ! -f "${INPUT_FILE}" ]]; then
        fail "File '${INPUT_FILE}' not found."
    fi
}

function usage()
{
    printf '\n\e[96m Helper script to create metadata files from unique database.\n\n Usage: '
    printf "\e[93m${BASH_SOURCE[0]} <option>\e[0m\n\n"
    printf '\e[96m Possible options:\n\n'
    printf '    \e[93m%-15s\e[0m  ->  \e[96m%s\e[0m\n' \
           '-a | --authors' 'Create AUTHORS.md metadata file'\
           '-c | --citation' 'Create CITATION.cff metadata file'\
           '-z | --zenodo'  'Create Zenodo .zenodo.json metadata file'
    printf '\n'
}

function regex_in_array()
{
    local element
    for element in "${@:2}"; do [[ "$element" =~ $1 ]] && return 0; done
    return 1
}

function parse_command_line_arguments()
{
    if regex_in_array '^-(h|-help)$' "$@"; then
        usage
        exit 0
    fi
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -a | --authors)
                CREATE_AUTHORS_FILE='TRUE'
                ;;
            -c | --citation)
                CREATE_CITATION_FILE='TRUE'
                ;;
            -z | --zenodo)
                CREATE_ZENODO_FILE='TRUE'
                ;;
            *)
                fail "Unrecognized option '${$1}'."
                ;;
        esac
        shift
    done
    if [[ ${CREATE_AUTHORS_FILE} = 'FALSE' && \
          ${CREATE_ZENODO_FILE} = 'FALSE' && \
          ${CREATE_CITATION_FILE} = 'FALSE' ]]; then
        warn "No output file was asked to be created, specify -a and/or -c and/or -z option(s)."
    fi
}

function create_authors_file_if_requested()
{
    if [[ ${CREATE_AUTHORS_FILE} = 'TRUE' ]]; then
        printf '\n Creating AUTHORS.md file...'
        {
            print_authors_header 'SMASH team'
            print_authors_section 'SMASH team'
            print_authors_header 'External contributors'
            print_authors_section 'External contributors'
            print_authors_header 'Acknowledgments'
            print_authors_section 'Acknowledgments'
        } > "${OUTPUT_AUTHORS_FILE}"
        printf ' done!\n'
    fi
}

function print_authors_header()
{
    local -r label="$1"
    case "${label}" in
        "SMASH team")
            cat <<EOF
# Authors

## The SMASH Team

Author  |  &ensp;E-Mail&ensp; | Copyright ©
 :----:  |  :----: | :---------
EOF
            ;;
        "External contributors")
            cat <<EOF

## External contributions

Year  | Author | E-Mail
:---: | :----: | :----:
EOF
            ;;
        Acknowledgments)
            cat <<EOF

## Further acknowledgments

Contributions by the following people are acknowledged.

Year  | Author | E-Mail
:---: | :----: | :----:
EOF
            ;;
        *)
            fail "Unrecognized label '${label}' in '${FUNCNAME[0]}' function."
            ;;
    esac
}

function print_authors_section()
{
    local -r label="$1"
    local type index number
    case "${label}" in
        'SMASH team')
            type='internal'
            ;;
        *)
            type='external'
            ;;
    esac
    number=$(jq ".\"${label}\" | length" "${INPUT_FILE}")
    for ((index=0; index<number; index++)); do
        print_${type}_author "$(jq ".\"${label}\"[${index}]" "${INPUT_FILE}")"
    done
}

function print_internal_author()
{
    local full_name email years
    set_author_variables "$1"
    printf '%s | %s | %s\n'\
        "${full_name}"\
        "[✉️](mailto:${email})"\
        "\`${years}\`"
}

function print_external_author()
{
    local full_name email years
    set_author_variables "$1"
    printf '%s | %s | %s\n'\
        "\`${years}\`"\
        "${full_name}"\
        "[✉️](mailto:${email})"
}

function set_author_variables()
{
    local -r json_author="$1"
    local first_name last_name
    first_name=$(jq -r '."first name"' <<< "${json_author}")
    last_name=$(jq -r '."last name"' <<< "${json_author}")
    full_name="${first_name} ${last_name}"
    email=$(jq -r '."email"' <<< "${json_author}")
    years=$(jq -r '."years"' <<< "${json_author}")
}


function create_zenodo_file_if_requested()
{
    if [[ ${CREATE_ZENODO_FILE} = 'TRUE' ]]; then
        printf '\n Creating .zenodo.json file...'
        # This command takes two JSON "files" and merge them letting possible identical
        # keys take the value of the last file (here no identical key should exist).
        # To make this work it is important to use the -s option of jq which reads all
        # inputs into an array and use it as the single input value.
        jq -s '.[0] + .[1]' \
            <(print_smash_team_for_zenodo) \
            <(print_related_identifiers_json_map) > "${OUTPUT_ZENODO_FILE}"
        printf ' done!\n'
    fi
}

function print_smash_team_for_zenodo()
{
    # This is a jq construct to basically get a different JSON file
    # from an existing one. It should not be difficult to understand it
    # referring to the well written jq documentation:
    #    https://jqlang.github.io/jq/manual/#basic-filters
    #
    # Some guidance:
    #   - The | operator is combining two filters (a kind of piping)
    #   - .XXX[] prints all entries of the JSON array named XXX
    #   - Quotes are needed when XXX contains spaces
    #   - The , operator concatenates two filters output which get the same input
    #   - The external {} and [] after creators: do an object and an array construction
    #   - The () group whatever is inside and around a key like (.XXX) make jq evaluate the key
    #   - The + for strings is a concatenation
    #   - The affiliation and orcid fields expand to the corresponding key:value pair
    #   - The sort_by and reverse do what the name says
    jq '
    {
        creators: [
            ."SMASH team"[], ."External contributors"[] |
            {name: (."last name"+", "+."first name"), affiliation, orcid}
        ] | sort_by(.name) | reverse
    }
    ' "${INPUT_FILE}"
}

function print_related_identifiers_json_map()
{
    cat <<EOF
{
  "related_identifiers": [
    {
      "scheme": "doi",
      "identifier": "https://doi.org/10.1103/PhysRevC.94.054905",
      "relation": "isDescribedBy",
      "resource_type": "publication-article"
    }
  ]
}
EOF
}

function create_citation_file_if_requested()
{
    if [[ ${CREATE_CITATION_FILE} = 'TRUE' ]]; then
        printf '\n Creating CITATION.cff file...'
        cat \
            <(print_citation_file_metadata) \
            <(print_authors_for_citation_file) \
            <(print_preferred_citation_for_citation_file) > "${OUTPUT_CITATION_FILE}"
        printf ' done!\n'
    fi
}

function print_citation_file_metadata()
{
    printf '%s\n' \
        'cff-version: 1.2.0' \
        'message: "If you use this software, please cite both the article from preferred-citation and the software itself."' \
        'title: "SMASH"' \
        "version: ${SMASH_VERSION}" \
        "date-released: ${SMASH_RELEASE_DATE}" \
        'identifiers:' \
        '  - type: doi' \
        '    value: 10.5281/zenodo.3484711' \
        '    description: "The concept DOI of the software."' \
        'license: GPL-3.0-or-later' \
        'repository-code: "https://github.com/smash-transport/smash"' \
        'url: "https://smash-transport.github.io"' \
        'type: software' \
        'contact:' \
        '  - email: elfner@itp.uni-frankfurt.de' \
        '    family-names: Elfner' \
        '    given-names: Hannah'\
        '    website: "https://www.elfner-group.science/index.html"' \
        ''
}

function print_preferred_citation_for_citation_file()
{
    printf '%s\n' \
        '' \
        'preferred-citation:' \
        '  title: "Particle production and equilibrium properties within a new hadron transport approach for heavy-ion collisions"' \
        '  authors:' \
        '    - family-names: Weil' \
        '      given-names: Janus' \
        '    - family-names: Steinberg' \
        '      given-names: Vinzent' \
        '    - family-names: Staudenmaier' \
        '      given-names: Jan' \
        '    - family-names: Pang' \
        '      given-names: Long-Gang' \
        '    - family-names: Oliinychenko' \
        '      given-names: Dmytro' \
        '    - family-names: Mohs' \
        '      given-names: Justin' \
        '    - family-names: Kretz' \
        '      given-names: Matthias' \
        '    - family-names: Kehrenberg' \
        '      given-names: Thomas' \
        '    - family-names: Goldschmidt' \
        '      given-names: Andy' \
        '    - family-names: Bäuchle' \
        '      given-names: Bjørn' \
        '    - family-names: Auvinen' \
        '      given-names: Jussi' \
        '    - family-names: Attems' \
        '      given-names: Maximilian' \
        '    - family-names: Petersen' \
        '      given-names: Hannah' \
        '  volume: 94' \
        '  issn: "2469-9993"' \
        '  url: "http://dx.doi.org/10.1103/PhysRevC.94.054905"' \
        '  doi: 10.1103/physrevc.94.054905' \
        '  number: 5' \
        '  journal: "Physical Review C"' \
        '  publisher: "American Physical Society (APS)"' \
        '  year: 2016' \
        '  month: 11' \
        '  type: article'
}

function print_authors_for_citation_file()
{
    # Refer to the comment in the 'print_smash_team_for_zenodo' function to get a
    # description of the following jq command. Here yq is used to convert the
    # JSON output to YAML format.
    jq '
    {
        authors: [
            ."SMASH team"[], ."External contributors"[] |
            {"given-names": ."first name", "family-names": ."last name", affiliation, orcid}] |
            sort_by(."family-names") | reverse
    }
    ' "${INPUT_FILE}" | yq -P -oy
}

function fail()
{
    printf "\n \e[1;91mERROR:\e[22m $1\e[0m\n" 1>&2
    shift
    while [[ $# -gt 0 ]]; do
        printf " \e[91m       $1\e[0m\n" 1>&2
        shift
    done
    exit 1
}

function warn()
{
    printf "\n \e[1;93mWARNING:\e[22m $*\e[0m\n"
    exit 1
}

#===================================================================

main "$@"
