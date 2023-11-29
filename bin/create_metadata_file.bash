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
# same place. The jq tool is used to deal with JSON files.
#===================================================================

# Bash stricter mode
set -euo pipefail
shopt -s inherit_errexit

readonly \
    INPUT_FILE='.authors.json'\
    OUTPUT_AUTHORS_FILE='AUTHORS.md'\
    OUTPUT_ZENODO_FILE='.zenodo.json'

CREATE_ZENODO_FILE='FALSE'
CREATE_AUTHORS_FILE='FALSE'

function main()
{
    make_preliminary_checks
    parse_command_line_arguments "$@"
    create_authors_file_if_requested
    create_zenodo_file_if_requested
}

function make_preliminary_checks()
{
    if ! hash jq; then
        fail\
            "Program 'jq' not found."\
            'Follow instructions at https://jqlang.github.io/jq/download/ to install it.'
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
           '-z | --zenodo'  'Create Zenodo .zenodo.json metadata file' \
           '-a | --authors' 'Create AUTHORS.md metadata file'
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
            -z | --zenodo)
                CREATE_ZENODO_FILE='TRUE'
                ;;
            *)
                fail "Unrecognized option '${$1}'."
                ;;
        esac
        shift
    done
    if [[ ${CREATE_AUTHORS_FILE} = 'FALSE' && ${CREATE_ZENODO_FILE} = 'FALSE' ]]; then
        warn "No output file was asked to be created, specify -a and/or -z option(s)."
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
