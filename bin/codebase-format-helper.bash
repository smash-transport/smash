#!/usr/bin/env bash

#===================================================
#
#    Copyright (c) 2018-2019,2021-2022
#      SMASH Team
#
#    GNU General Public License (GPLv3 or later)
#
#===================================================

shopt -s globstar
CHOSEN_LANGUAGE=''
FORMATTER_COMMAND=''
MODE_PERFORM='FALSE'
MODE_TEST='FALSE'
FILES_TO_FORMAT=()

function main()
{
    parse_command_line_arguments "$@"
    check_formatter_availability
    check_formatter_version
    look_for_files_to_format
    if [[ ${MODE_PERFORM} = 'TRUE' ]]; then
        perform_formatting
    elif [[ ${MODE_TEST} = 'TRUE' ]]; then
        test_formatting
    fi
}

#===================================================================
# Auciliary functions

function check_formatter_availability()
{
    if ! type "${FORMATTER_COMMAND}" &>/dev/null; then
        fail "${FORMATTER_COMMAND} command not found."
    fi
}

function check_formatter_version()
{
    local required found
    found=$(${FORMATTER_COMMAND} --version)
    if [[ ${CHOSEN_LANGUAGE} = 'C++' ]]; then
        required='6.0.0'
        found="${found:21:5}"
    elif [[ ${CHOSEN_LANGUAGE} = 'CMake' ]]; then
        required='0.6.13'
    fi
    if [[ "${found}" != "${required}" ]]; then
        fail "Wrong ${FORMATTER_COMMAND} version found: ${found} (${required} is required)."
    fi
}

function look_for_files_to_format()
{
    local base_dir
    base_dir="$(dirname $BASH_SOURCE[0])/.."
    if [[ ${CHOSEN_LANGUAGE} = 'C++' ]]; then
        FILES_TO_FORMAT=(
            "${base_dir}/src"/**/*.{cc,h}
        )
    elif [[ ${CHOSEN_LANGUAGE} = 'CMake' ]]; then
        FILES_TO_FORMAT=(
            "${base_dir}"/**/CMakeLists.txt
            "${base_dir}"/**/*.cmake
        )
        for index in "${!FILES_TO_FORMAT[@]}"; do
            if [[ ${FILES_TO_FORMAT[index]} =~ ^${base_dir}/(3rdparty|build[^/]*)/ ]]; then
                unset -v 'FILES_TO_FORMAT[index]'
            fi
        done
        FILES_TO_FORMAT=("${FILES_TO_FORMAT[@]}")
    fi
}

function perform_formatting()
{
    local filename
    for filename in "${FILES_TO_FORMAT[@]}"; do
        printf "Formatting ${filename}\n"
        ${FORMATTER_COMMAND} -i "${filename}"
    done
}

function test_formatting()
{
    local check_flag file_under_examination formatting_differences
    printf "Testing that ${FORMATTER_COMMAND} does not change the source code in the working directory...\n"
    check_flag=0
    for file_under_examination in "${FILES_TO_FORMAT[@]}"; do
        # NOTE: cmake-format has a --check option implemented that might be used
        #       here, but if there are differences we want to print them to the
        #       user and hence we just store them here avoiding to format twice
        formatting_differences=$(diff <(cat "${file_under_examination}") <(${FORMATTER_COMMAND} "${file_under_examination}") 2>&1)
        if [[ $? -ne 0 ]]; then
            printf "File \"${file_under_examination}\" not properly formatted. Comparison with a properly formatted file:\n"
            printf "${formatting_differences}\n"
            check_flag=1
        fi
    done
    if [[ ${check_flag} -eq 1 ]]; then
        fail "${FORMATTER_COMMAND} was not properly run on latest commit."
    else
        printf "\n\e[1;92m PASS:\e[22m No changes to source code by ${FORMATTER_COMMAND}.\e[0m\n\n"
    fi
}

function usage()
{
    printf '\n\e[96m'
    printf '%s\n' \
           ' Helper script to format source files in the codebase.' \
           '' \
           " Usage: ${BASH_SOURCE[0]} (C++|CMake) <option>" \
           '' \
           ' Possible options:'
    printf '    %-15s  ->  %s\n' \
           '-p | --perform' 'Perform automatic formatting (for developers)' \
           '-t | --test' 'Test that automatic formatting has been performed (for CI builds)' \
           '-h | --help' 'Display this hel'
    printf '\n\e[0m'
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
    if [[ ! $1 =~ ^C(\+\+|Make)$ ]]; then
        fail "First command line option must be either 'C++' or 'CMake'."
    else
        CHOSEN_LANGUAGE=$1
        shift
    fi
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -p | --perform)
                MODE_PERFORM='TRUE'
                ;;
            -t | --test)
                MODE_TEST='TRUE'
                ;;
            *)
                fail "Unrecognised option."
                ;;
        esac
        shift
    done
}

function fail()
{
    printf "\n \e[1;91mERROR:\e[22m $*\e[0m\n\n" 1>&2
    exit 1
}

#===================================================================

main "$@"
