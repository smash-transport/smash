#!/usr/bin/env bash

#===================================================
#
#    Copyright (c) 2022-2024
#      SMASH Team
#
#    GNU General Public License (GPLv3 or later)
#
#===================================================

# NOTE: We use the 'BM_' prefix for global variables and in particular the idea
#       behind this script is to extract the benchmark timings with their errors
#       from the benchmark files into global associative arrays which have the
#       names of the benchmark as keys.

trap 'printf "\n"' EXIT

function main()
{
    printf '\n'
    Validate_command_line_options "$@"
    Extract_Benchmarks_Timings_Into_Global_Arrays "$@"
    if [[ $# -eq 2 ]]; then
        Compare_Two_Benchmarks "$1" "$2"
    fi
}

function Validate_command_line_options()
{
    if [[ $# -lt 2 ]]; then
        printf 'Usage: %s file_1 file_2 [further_files]...\n' "${BASH_SOURCE}"
        exit 1
    fi
    local filename
    for filename in "$@"; do
        if [[ ! -f "${filename}" ]]; then
            printf "File '${filename}' not found.\n"
            exit 1
        elif [[ $(basename "${filename}") != bm-results-SMASH-*.md ]]; then
            printf '%s\n'\
                   "File '${filename}' does not seem to be a benchmark result file."\
                   "Benchmark results files must match the 'bm-results-SMASH-*.md' glob."
            exit 1
        fi
    done
}

function Extract_Benchmarks_Timings_Into_Global_Arrays()
{
    local filename counter=1
    for filename in "$@"; do
        declare -g -A BM_timings_${counter}{,_errors}
        Extract_benchmarks_data_from_file "${filename}" into_array BM_timings_${counter}
        (( counter++ ))
    done
}

function Compare_Two_Benchmarks()
{
    readonly version_old=$(Extract_version_from_filename "$1")
    readonly version_new=$(Extract_version_from_filename "$2")
    # Print report
    local old_dismissed_benchmarks new_introduced_benchmarks benchmark
    old_dismissed_benchmarks=()
    new_introduced_benchmarks=()
    Set_length_longest_label "${!BM_timings_1[@]}" # TODO: Adjust if longest was dismissed
    Print_report_header
    for benchmark in "${!BM_timings_1[@]}"; do
        if [[ ${BM_timings_2["${benchmark}"]} = '' ]]; then
            old_dismissed_benchmarks+=( "${benchmark}" )
        else
            Print_report_line "${benchmark}" ${BM_timings_1["${benchmark}"]} ${BM_timings_2["${benchmark}"]}
        fi
    done
    for benchmark in "${!BM_timings_2[@]}"; do
        if [[ ${BM_timings_1["${benchmark}"]} = '' ]]; then
            new_introduced_benchmarks+=( "${benchmark}" )
        fi
    done
    Print_list_of_benchmarks "dismissed after ${version_old}" "${old_dismissed_benchmarks[@]}"
    Print_list_of_benchmarks "newly introduced ${version_new}" "${new_introduced_benchmarks[@]}"
}

#========================== Utility functions ==================================

function Extract_version_from_filename()
{
    local filename="${1##*bm-results-}"
    printf "${filename%.md}"
}

# Function requires to be invoked with 'into_array' as second argument
function Extract_benchmarks_data_from_file()
{
    local filename array_name errors_array_name
    filename=$1
    array_name=$3
    errors_array_name="${array_name}_errors"
    if [[ $2 != 'into_array' ]]; then
        printf "Function '${FUNCNAME}' wrongly called.\n"
        exit 2
    elif [[ $3 =~ [^a-zA-Z0-9_] ]]; then
        printf "Function '${FUNCNAME}' must be called with valid array_name as third argument.\n"
        exit 1
    elif ! declare -p $3 &> /dev/null; then
        printf "Array '$3' must be declared before calling '${FUNCNAME}'.\n"
        exit 1
    elif ! declare -p "${3}_errors" &> /dev/null; then
        printf "Array '${3}_errors' must be declared before calling '${FUNCNAME}'.\n"
        exit 1
    fi
    local list_of_benchmarks list_of_timings list_of_errors
    readarray -t list_of_benchmarks < <(grep -E '^### [A-Z]' "${filename}" | sed 's/### //')
    readarray -t list_of_timings    < <(grep -E 'seconds time elapsed' "${filename}" | awk '{print $1}')
    readarray -t list_of_errors     < <(grep -E 'seconds time elapsed' "${filename}" | awk '{print $3}')
    for index in "${!list_of_benchmarks[@]}"; do
        eval ${array_name}["'"${list_of_benchmarks[index]}"'"]=${list_of_timings[index]}
        eval ${errors_array_name}["'"${list_of_benchmarks[index]}"'"]=${list_of_errors[index]}
    done
}

function Set_length_longest_label()
{
    local longest_label_length=0
    for label in "$@"; do
        if [[ ${#label} -gt ${longest_label_length} ]]; then
            longest_label_length=${#label}
        fi
    done
    readonly width_benchmark_column="${longest_label_length}"
}

function Print_report_header()
{
    printf "%${width_benchmark_column}s%25s%25s%18s\n"\
           'BENCHMARK' "${version_old} {s}" "${version_new} {s}" 'Time change'
}

function Print_report_line()
{
    local delta color
    delta="$(awk -v old="$2" -v new="$3" 'BEGIN{printf "%+6.2f", (new-old)/new*100}')"
    color="$(awk -v x="${delta}" 'BEGIN{printf "%d", (x>0) ? 91 : 92}')"
    printf "\e[96m%${width_benchmark_column}s\e[0m%25s%25s\e[${color}m%18s\e[0m\n" "$1" "$2" "$3" "${delta}%"
}

function Print_list_of_benchmarks()
{
    local label=$1
    shift
    if [[ $# -gt 0 ]]; then
        printf '\nThe following benchmarks have been %s\n' "${label}"
        for benchmark in "$@"; do
            printf '  - %s\n' "${benchmark}"
        done
    fi
}


main "$@"
