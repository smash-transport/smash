#!/usr/bin/env bash

trap 'printf "\n"' EXIT

function main()
{
    printf '\n'
    Validate_command_line_options "$@"
    readonly results_filename_old="$1"
    readonly results_filename_new="$2"
    readonly version_old=$(Extract_version_from_filename "${results_filename_old}")
    readonly version_new=$(Extract_version_from_filename "${results_filename_new}")
    # Parse benchmark results
    declare -A benchmarks_old benchmarks_new
    Extract_benchmarks_data_from_file "${results_filename_old}" into_array 'benchmarks_old'
    Extract_benchmarks_data_from_file "${results_filename_new}" into_array 'benchmarks_new'
    # Print report
    old_dismissed_benchmarks=()
    new_introduced_benchmarks=()
    Set_length_longest_label "${!benchmarks_old[@]}" # TODO: Adjust if longest was dismissed
    Print_report_header
    for benchmark in "${!benchmarks_old[@]}"; do
        if [[ ${benchmarks_new["${benchmark}"]} = '' ]]; then
            old_dismissed_benchmarks+=( "${benchmark}" )
        else
            Print_report_line "${benchmark}" ${benchmarks_old["${benchmark}"]} ${benchmarks_new["${benchmark}"]}
        fi
    done
    for benchmark in "${!benchmarks_new[@]}"; do
        if [[ ${benchmarks_old["${benchmark}"]} = '' ]]; then
            new_introduced_benchmarks+=( "${benchmark}" )
        fi
    done
    Print_list_of_benchmarks "dismissed after ${version_old}" "${old_dismissed_benchmarks[@]}"
    Print_list_of_benchmarks "newly introduced ${version_new}" "${new_introduced_benchmarks[@]}"
}
 
function Validate_command_line_options()
{
    if [[ $# -ne 2 ]]; then
        printf 'Usage: %s <old_benchmark_results_file> <new_benchmark_results_file>\n' "${BASH_SOURCE}"
        exit 1
    fi
    if [[ ! -f $1 ]]; then
        printf "File '$1' not found.\n"
        exit 1
    elif [[ $(basename "$1") != bm-results-SMASH-*.md ]]; then
        printf "File '$1' does not seem to be a benchmark result file.\n"
        exit 1
    fi
    if [[ ! -f $2 ]]; then
        printf "File '$2' not found.\n"
        exit 1
    elif [[ $(basename "$2") != bm-results-SMASH-*.md ]]; then
        printf "File '$1' does not seem to be a benchmark result file.\n"
        exit 1
    fi
}

function Extract_version_from_filename()
{
    local filename="${1##*bm-results-}"
    printf "${filename%.md}"
}

# Function requires to be invoked with 'into_array' as second argument
function Extract_benchmarks_data_from_file()
{
    if [[ $2 != 'into_array' ]]; then
        printf "Function '${FUNCNAME}' wrongly called.\n"
        exit 2
    elif [[ $3 =~ [^a-zA-Z_] ]]; then
        printf "Function '${FUNCNAME}' must be called with valid array_name as third argument.\n"
        exit 1        
    fi
    local filename array_name list_of_benchmarks list_of_timings
    filename=$1
    array_name=$3
    readarray -t list_of_benchmarks < <(grep -E '^### [A-Z]' "${filename}" | sed 's/### //')
    readarray -t list_of_timings    < <(grep -E 'seconds time elapsed' "${filename}" | awk '{print $1}')
    for index in "${!list_of_benchmarks[@]}"; do
        eval ${array_name}["'"${list_of_benchmarks[index]}"'"]=${list_of_timings[index]}
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
    printf "%${width_benchmark_column}s%20s%20s%15s\n" 'BENCHMARK' "${version_old} {s}" "${version_new} {s}" 'Time change'
}

function Print_report_line()
{
    local delta color
    delta="$(awk -v old="$2" -v new="$3" 'BEGIN{printf "%+6.2f", (new-old)/new*100}')"
    color="$(awk -v x="${delta}" 'BEGIN{printf "%d", (x>0) ? 91 : 92}')"
    printf "\e[96m%${width_benchmark_column}s\e[0m%20s%20s\e[${color}m%15s\e[0m\n" "$1" "$2" "$3" "${delta}%"
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
