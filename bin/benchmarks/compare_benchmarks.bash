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
    Parse_command_line_options "$@"
    printf '\n'
    Extract_Benchmarks_Timings_Into_Global_Arrays "${BM_benchmarks_files[@]}"
    if [[ $# -eq 2 ]]; then
        Compare_Two_Benchmarks "${BM_benchmarks_files[@]}"
    fi
    Make_Benchmark_Plot "${BM_benchmarks_files[@]}"
}

function Parse_command_line_options()
{
    local option filename
    for option in "$@"; do
        if [[ ${option} =~ ^-(h|-help)$ ]]; then
            printf '\nUsage: \e[96m%s [OPTIONS] [--] file_1 file_2 [further_files]...\e[0m\n\nAccepted options:\n' "${BASH_SOURCE}"
            printf '\e[93m%20s\e[0m  ->  \e[96m%s\e[0m\n'\
                '-a | --annotate' 'Which bars to annotate as Python bool or list (default: True).'
            exit 0
        fi
    done
    BM_bars_annotations='True'
    BM_benchmarks_files=()
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -a | --annotate )
                if [[ $2 =~ ^(True|False|\[[0-9]+(,[0-9]+)*\])$ ]]; then
                    BM_bars_annotations="$2"
                else
                    Fail "Error specifying value of '-a' option.\n"
                fi
                shift 2
                ;;
            -- )
                shift
                BM_benchmarks_files+=( "$@" )
                shift $#
                ;;
            [^-]* )
                BM_benchmarks_files+=( "$1" )
                shift 1
                ;;
            * )
                Fail "Unrecognized option '$1'."
                ;;
        esac
    done
    for filename in "${BM_benchmarks_files[@]}"; do
        if [[ ! -f "${filename}" ]]; then
            Fail "File '${filename}' not found."
            exit 1
        elif [[ $(basename "${filename}") != bm-results-SMASH-*.md ]]; then
            Fail \
                "File '${filename}' does not seem to be a benchmark result file."\
                "Benchmark results files must match the 'bm-results-SMASH-*.md' glob."
            exit 1
        fi
    done
    if [[ "${#BM_benchmarks_files[@]}" -lt 2 ]]; then
        Fail 'At least two benchmark files must be passed to the script.'
    fi
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

# Here we use the python plotting script in the same folder for which data
# has to be prepared in a very precise form -> Do it in constructing the
# python dictionary as bash string and then inject it in the python code
#
# NOTE: Since it is in general not guaranteed that all benchmark files contain
#       the same SMASH runs, here first all existing runs are extracted in
#       a list and then we go through such a list to put data in the correct
#       format (cf. with the python plotting script).
#
# NOTE: The data are extracted into global associative arrays and here it is
#       needed to iterate over them, which is not a native feature of bash.
#       Hence we create a reference to each array in a loop and use it as needed.
function Make_Benchmark_Plot()
{
    local data_dictionary x_ticks number_of_arrays list_of_benchmarks index name
    number_of_arrays=$#
    list_of_benchmarks=()
    for ((index=1; index<=number_of_arrays; index++)); do
        for name in BM_timings_${index}{,_errors}; do
            declare -n array_ref="${name}"
            list_of_benchmarks+=( "${!array_ref[@]}" )
        done
    done
    # Remove duplicates from the list and then iterate through it
    readarray -d $'\0' -t list_of_benchmarks < <(printf "%s\0" "${list_of_benchmarks[@]}" | sort -uz)
    data_dictionary=''
    for ((index=1; index<=number_of_arrays; index++)); do
        data_dictionary+="  \"$(Extract_version_from_filename "${!index}")\": ["
        for name in "${list_of_benchmarks[@]}"; do
            declare -n \
                value_array_ref="BM_timings_${index}"\
                errors_array_ref="BM_timings_${index}_errors"
            if [[ ${data_dictionary: -1} != '[' ]]; then
                data_dictionary+=', '
            fi
            if [[ ${value_array_ref["${name}"]} != '' ]]; then
                data_dictionary+="(${value_array_ref["${name}"]},${errors_array_ref["${name}"]})"
            else
                data_dictionary+="(0,0)"
            fi
        done
        data_dictionary+=$'],\n'
    done
    data_dictionary='data = {'$'\n'"${data_dictionary%,$'\n'}"$'\n''}'
    printf -v x_ticks '"%s", ' "${list_of_benchmarks[@]}"
    x_ticks="[${x_ticks%, }]"
    # Shorten names for x-ticks
    x_ticks="${x_ticks// Run/}"             # remove ' Run'
    x_ticks="${x_ticks//without/w\/o}"      # 'without' -> 'w/o'
    x_ticks="${x_ticks//with/w.}"           # 'with' -> 'w.'
    x_ticks="${x_ticks//Testparticles/TP}"  # 'Testparticles' -> 'TP'
    # Make final plot
    python3 - <<EOF
from multi_bar_histogram_plot import make_multi_bar_histogram_plot

${data_dictionary}

make_multi_bar_histogram_plot(data, total_width=.75, single_width=.95,
                              xticks=${x_ticks},
                              ylabel="Run time {s}",
                              ymargin=0.2,
                              annotate=${BM_bars_annotations})
EOF
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

function Fail()
{
    printf "\n \e[1;91mERROR:\e[22m $*\e[0m\n" 1>&2
    exit 1
}

main "$@"
