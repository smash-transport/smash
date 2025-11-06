#!/usr/bin/env bash

#=================================================
#
#    Copyright (c) 2023,2025
#      SMASH Team
#
#    GNU General Public License (GPLv3 or later)
#
#=================================================

trap 'printf "\n"' EXIT

shopt -s extglob

declare -rg minimum_bash_version='4.3.0'

smash_sources_root_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )/.." &> /dev/null && pwd )

# Current old Pythia version (automatically detected)
current_pythia_version=''

# New Pythia version (the script asks the user)
new_pythia_version=''

#=================================================
# Principal functions (i.e. main and the functions directly called by main)

function main()
{
    check_bash_version
    find_current_pythia_version
    ask_new_pythia_version
    ask_confirmation
    change_pythia_version
    print_closing_message
    # This function could also be called by another function, but if we arrive here we're done and we can exit
    exit $?
}

function check_bash_version()
{
    local required found
    required=${minimum_bash_version}
    found=$(echo ${BASH_VERSINFO[@]:0:3} | tr ' ' '.')
    if [[ $(printf '%s\n' "${required}" "${found}" | sort -V | head -n1) != "${required}" ]]; then
        printf "Minimum bash version required is ${required}, but found ${found}.\n"
        inform_about_ending_and_exit
    fi
}

function find_current_pythia_version()
{
    local pythia_version_in_readme pythia_version_in_cmakelist
    pythia_version_in_readme=$(grep 'Pythia' ${smash_sources_root_dir}/README.md | grep -Eo '[1-9]?[0-9]\.[0-9]{3}'\
        | sort -u)
    pythia_version_in_cmakelist=$(grep 'Pythia' ${smash_sources_root_dir}/src/CMakeLists.txt | grep -Eo '[1-9]?[0-9]\.[0-9]{3}'\
        | sort -u)
    if [[ "${pythia_version_in_readme}" != "${pythia_version_in_cmakelist}" ]]
    then
        printf "Different Pythia versions detected in README.md (%s) and in src/CMakeLists.txt (%s).\n"\
               ${pythia_version_in_readme} ${pythia_version_in_cmakelist}
        printf "The automatic Pythia version update with this script is aborted because it would be unsafe.\n"
        inform_about_ending_and_exit
    else
        printf "Detected Pythia version ${pythia_version_in_readme}.\n\n"
        current_pythia_version=${pythia_version_in_readme}
    fi
}

function ask_new_pythia_version()
{
    local tmp_new_pythia_version
    printf "Please, insert the new Pythia version to which you want to upgrade. (\e[1;36me \e[0mfor \e[1;36mexit\e[0m): "
    read tmp_new_pythia_version
    case "${tmp_new_pythia_version}" in
        "e")
            inform_about_ending_and_exit 0
            ;;
        [1-9]?([0-9]).[0-9][0-9][0-9])
            new_pythia_version=${tmp_new_pythia_version}
            ;;
        *)
            printf "\e[1;31mIncorrect (or missing) Pythia version number format\e[0m\n"
            printf "Expecting something like x.zkj or xy.zkj, where y,z,k and j are digits from 0 to 9 and x from 1 to 9\n\n"
            ask_new_pythia_version
            ;;
    esac
}

function ask_confirmation()
{
    printf "Replacing occurrences of Pythia version %s with version %s.\n" ${current_pythia_version} ${new_pythia_version}
    printf "Do you want to continue?\n"
    if user_said_no
    then
       printf "OK. We restart from the beginning.\n\n"
       main
    fi
}

function change_pythia_version()
{
    local current_pythia_version_no_dots new_pythia_version_no_dots files_to_be_processed \
          current_pythia_main_version new_pythia_main_version
    current_pythia_version_no_dots=${current_pythia_version/./}
    current_pythia_main_version=${current_pythia_version_no_dots::-2}
    new_pythia_version_no_dots=${new_pythia_version/./}
    new_pythia_main_version=${new_pythia_version_no_dots::-2}

    files_to_be_processed=(
                           "${smash_sources_root_dir}"/README.md
                           "${smash_sources_root_dir}"/INSTALL.md
                           "${smash_sources_root_dir}"/.github/workflows/**.yml
                           "${smash_sources_root_dir}"/cmake/FindSMASH.cmake
                           "${smash_sources_root_dir}"/containers/**/Dockerfile
                           "${smash_sources_root_dir}"/examples/**/README.md
                           "${smash_sources_root_dir}"/src/CMakeLists.txt
                           )
    # NOTE: We use -i'.bkp' here so that this script works with both GNU and BSD sed implementations.
    #       The only drawback is that we need to remove '.bkp' files if sed succeeds.
    sed -i'.bkp' -r -e 's/([pP]ythia.*)'"${current_pythia_version}"'/\1'"${new_pythia_version}"'/g'\
        -e 's/([pP]ythia.*)'"${current_pythia_version_no_dots}"'/\1'"${new_pythia_version_no_dots}"'/g'\
        -e 's/pythia'"${current_pythia_version_no_dots}"'/pythia'"${new_pythia_version_no_dots}"'/g'\
        -e 's/download\/pythia'"${current_pythia_main_version}"'/download\/pythia'"${new_pythia_main_version}"'/g'\
        "${files_to_be_processed[@]}"
    if [[ $? -eq 0 ]]; then
        rm "${files_to_be_processed[@]/%/.bkp}"
    else
        printf "\nThe sed command to automatically change Pythia version failed! Leaving backup files, please investigate!\n"
        inform_about_ending_and_exit
    fi
}

#=================================================
# Auxiliary functions (i.e. not directly called by the main function)

function inform_about_ending_and_exit()
{
    printf "\n \e[1;91mEXITING THE SCRIPT\e[0m\n" 1>&2
    exit ${1-1}
}

function user_said_yes()
{
    local user_answer
    printf " \e[1;32my (YES, default option)\e[0m OR \e[1;31mn (NO) \e[0m? (\e[1;36me \e[0mfor \e[1;36mexit\e[0m): "
    while read user_answer; do
        case ${user_answer} in
            "y")
                return 0
                ;;
            "n")
                return 1
                ;;
            "e")
                inform_about_ending_and_exit 0
                ;;
            *)
                printf "\e[1;95mInvalid choice!\e[0m\n\n"
                ;;
        esac
    done
}

function user_said_no()
{
    ! user_said_yes
}

function print_closing_message()
{
    printf "\n\e[35mPythia version changed. \e[1;35mPlease, remember to mention this operation in the file CHANGELOG.md!\e[0m\n"
}
#=================================================
# Main function
main
