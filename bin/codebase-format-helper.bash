#!/usr/bin/env bash

#===================================================
#
#    Copyright (c) 2018-2019,2021-2023
#      SMASH Team
#
#    GNU General Public License (GPLv3 or later)
#
#===================================================

trap 'printf "\n"' EXIT

#===================================================================
# Required versions of formatting programs (global variables)

declare -rg minimum_bash_version='4.3.0'
declare -rg clang_format_required_version='13.0'
declare -rg cmake_format_required_version='0.6.13'

#===================================================================
# Principal functions

function main()
{
    check_bash_version
    do_variables_and_shell_options_setup
    parse_command_line_arguments "$@"
    format_chosen_languages
}

function do_variables_and_shell_options_setup()
{
    shopt -s globstar nullglob
    CHOSEN_LANGUAGES=''
    declare -gA FORMATTER_COMMAND=(
        ['C++']='clang-format-'${clang_format_required_version%%.*}
        ['CMake']='cmake-format'
    )
    MODE_PERFORM='FALSE'
    MODE_TEST='FALSE'
    VERBOSE='TRUE'
    FILES_TO_FORMAT=()
}

function check_bash_version()
{
    local required found
    required=${minimum_bash_version}
    found=$(echo ${BASH_VERSINFO[@]:0:3} | tr ' ' '.')
    if [[ $(printf '%s\n' "${required}" "${found}" | sort -V | head -n1) != "${required}" ]]; then
        fail "Minimum bash version required is ${required}, but found ${found}."
    fi
}

function format_chosen_languages()
{
    local language
    for language in "${CHOSEN_LANGUAGES[@]}"; do
        format_given_language "${language}"
    done
}

function format_given_language()
{
    check_formatter_availability "$1"
    check_formatter_version "$1"
    look_for_files_to_format "$1"
    perform_or_test_formatting "$1"
}

function perform_or_test_formatting()
{
    if [[ ${MODE_PERFORM} = 'TRUE' ]]; then
        perform_formatting "$1"
    elif [[ ${MODE_TEST} = 'TRUE' ]]; then
        test_formatting "$1"
    fi
}

#===================================================================
# Auxiliary functions

function check_formatter_availability()
{
    if ! hash "${FORMATTER_COMMAND[$1]}" &>/dev/null; then
        if [[ $1 = 'C++' ]] && [[ "${FORMATTER_COMMAND[$1]}" != 'clang-format' ]]; then
            FORMATTER_COMMAND['C++']='clang-format'
            check_formatter_availability 'C++'
        else
            fail "'${FORMATTER_COMMAND[$1]}' command not found."
        fi
    else
        readonly -A FORMATTER_COMMAND[$1]
    fi
}

function check_formatter_version()
{
    local language required found
    language="$1"
    found=$(${FORMATTER_COMMAND[${language}]} --version)
    if [[ ${language} = 'C++' ]]; then
        required=${clang_format_required_version}
        found=$(grep -oE '[1-9][0-9]*[.][0-9]+' <<< "${found}" | head -n1)
    elif [[ ${language} = 'CMake' ]]; then
        required=${cmake_format_required_version}
    fi
    if [[ "${found}" != "${required}" ]]; then
        fail "Wrong ${FORMATTER_COMMAND[${language}]} version found: ${found} (${required} is required)."
    fi
}

function look_for_files_to_format()
{
    local language base_dir
    language="$1"
    base_dir="$(dirname ${BASH_SOURCE[0]})/.."
    if [[ ${language} = 'C++' ]]; then
        # C++ extenstions accepted by GNU compiler: https://gcc.gnu.org/onlinedocs/gcc/Overall-Options.html
        FILES_TO_FORMAT=(
            "${base_dir}/"{src,examples}/**/*.{h,hh,H,hp,hxx,hpp,HPP,h++,tcc,cc,cp,cxx,cpp,CPP,c++,C,ii}
        )
    elif [[ ${language} = 'CMake' ]]; then
        FILES_TO_FORMAT=(
            "${base_dir}"/**/CMakeLists.txt
            "${base_dir}"/**/*.cmake
        )
        for index in "${!FILES_TO_FORMAT[@]}"; do
            if [[ ${FILES_TO_FORMAT[index]} =~ ^${base_dir}/3rdparty/(Cuba[^/]*/|CMakeLists.txt$) ]]; then
                continue
            elif [[ ${FILES_TO_FORMAT[index]} =~ ^${base_dir}/(3rdparty|build[^/]*)/ ]]; then
                unset -v 'FILES_TO_FORMAT[index]'
            fi
        done
        FILES_TO_FORMAT=("${FILES_TO_FORMAT[@]}")
    fi
}

function perform_formatting()
{
    local language filename
    language="$1"
    verbose ''
    for filename in "${FILES_TO_FORMAT[@]}"; do
        verbose "Formatting ${filename}"
        ${FORMATTER_COMMAND[${language}]} -i "${filename}"
    done
    printf "\n\e[1;92m DONE:\e[22m ${language} code formatted by ${FORMATTER_COMMAND[${language}]}.\e[0m\n"
}

function test_formatting()
{
    local language check_flag file_under_examination formatting_differences
    language="$1"
    verbose "\n Testing that ${FORMATTER_COMMAND[${language}]} does not change the source code in the working directory..."
    check_flag=0
    for file_under_examination in "${FILES_TO_FORMAT[@]}"; do
        # NOTE: cmake-format has a --check option implemented that might be used
        #       here, but if there are differences we want to print them to the
        #       user and hence we just store them here avoiding to format twice
        formatting_differences=$(diff <(cat "${file_under_examination}") <(${FORMATTER_COMMAND[${language}]} "${file_under_examination}") 2>&1)
        if [[ $? -ne 0 ]]; then
            verbose "\n File \"${file_under_examination}\" not properly formatted. Comparison with a properly formatted file:"
            verbose "\e[0m${formatting_differences}"
            check_flag=1
        fi
    done
    if [[ ${check_flag} -eq 1 ]]; then
        fail "${FORMATTER_COMMAND[${language}]} was not properly run on latest commit."
    else
        printf "\n\e[1;92m PASS:\e[22m No changes to source code by ${FORMATTER_COMMAND[${language}]}.\e[0m\n"
    fi
}

function usage()
{
    printf '\n\e[96m Helper script to format source files in the codebase.\n\n Usage: '
    printf "\e[93m${BASH_SOURCE[0]} [C++|CMake] <option>\e[0m\n\n"
    printf '\e[96m Possible options:\n\n'
    printf '    \e[93m%-15s\e[0m  ->  \e[96m%s\e[0m\n' \
           '-p | --perform' 'Perform automatic formatting (for developers)' \
           '-t | --test' 'Test that automatic formatting has been performed (for CI builds)' \
           '-q | --quiet' 'Reduce output to minimum' \
           '-h | --help' 'Display this help'
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
    if [[ ! $1 =~ ^C(\+\+|Make)$ ]]; then
        CHOSEN_LANGUAGES=( 'C++' 'CMake' )
    else
        CHOSEN_LANGUAGES=( "$1" )
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
            -q | --quiet)
                VERBOSE='FALSE'
                ;;
            *)
                fail "Unrecognised option."
                ;;
        esac
        shift
    done
    if [[ ${MODE_PERFORM} = 'FALSE' && ${MODE_TEST} = 'FALSE' ]]; then
        fail 'Either -t or -p option should be given. Run with -h for more information.'
    fi
}

function fail()
{
    printf "\n \e[1;91mERROR:\e[22m $*\e[0m\n" 1>&2
    exit 1
}

function verbose()
{
    if [[ ${VERBOSE} = 'TRUE' ]]; then
        printf " \e[96m$*\e[0m\n"
    fi
}

#===================================================================

main "$@"
