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
# Global variables

# Minimum required bash version (constant)
declare -rg minimum_bash_version='4.3.0'

# SMASH source directory (the script asks the user, initially guessing ../)
declare -g smash_sources_root_dir="$(cd .. && echo $PWD)"

# Current old Pythia version (automatically detected)
declare -g current_pythia_version=''

# New Pythia version (the script asks the user)
declare -g new_pythia_version=''

#=============================================================================
# Principal functions (i.e. main and the functions directly called by main)

function main()
{
    check_bash_version
    confirm_smash_source_directory
    find_current_pythia_version
    ask_new_pythia_version
    ask_confirmation
    change_pythia_version
    print_closing_message
    # This function could also be called by another function, but if we arrive here we're done and we can exit
    exit 0
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

function confirm_smash_source_directory()
{
    local choice
    printf "Is \e[33m${smash_sources_root_dir}\e[0m the root directory of the SMASH sources that you want to update?\n"
    ask_yes_or_no
    choice=$? 
    if [ ${choice} -ne 0 ]
    then
        ask_user_smash_source_tree_path
    fi
}

function find_current_pythia_version()
{
    local pythia_version_in_readme pythia_version_in_cmakelist
    pythia_version_in_readme=$(grep Pythia ${smash_sources_root_dir}/README.md | grep -Po '[1-9]?[0-9]\.[0-9]{3}')
    pythia_version_in_cmakelist=$(grep Pythia ${smash_sources_root_dir}/src/CMakeLists.txt | grep -Po '[1-9]?[0-9]\.[0-9]{3}')
    if [ "${pythia_version_in_readme}" != "${pythia_version_in_cmakelist}" ]
    then
        printf "Different Pythia versions detected in README.md (%s) and in src/CMakeLists.txt (%s).\n"\
               ${pythia_version_in_readme} ${pythia_version_in_cmakelist}
        printf "The automatic Pythia version update with this script is aborted because it would be unsafe.\n"
        pretty_exit
    else
        printf "Detected Pythia version ${pythia_version_in_readme}.\n\n"
        current_pythia_version=${pythia_version_in_readme}
    fi
}

function ask_new_pythia_version()
{
    local tmp_new_pythia_version
    printf "Please, insert the new Pythia version to which you want to upgrade. (\e[1;36me \e[0mfor \e[1;36mexit\e[0m)\n"
    read tmp_new_pythia_version
    if [ "${tmp_new_pythia_version}" == "e" ]
    then
    pretty_exit
    fi
    if ! [[ "${tmp_new_pythia_version}" =~ [1-9]?[1-9]\.[0-9]{3} ]]
    then
        printf "Incorrect (or missing) Pythia version number format\n"
        printf "Expecting something like x.zkj or xy.zkj, where y,z,k and j are digits from 0 to 9 and x from 1 to 9\n"
        pretty_exit
    else
        new_pythia_version=${tmp_new_pythia_version}
    fi
}

function ask_confirmation()
{
    local choice
    printf "Replacing occurrences of Pythia version %s with version %s.\n" ${current_pythia_version} ${new_pythia_version}
    printf "Do you want to continue?\n"
    ask_yes_or_no
    choice=$?
    if [ ${choice} -ne 0 ]
    then
       printf "OK. We restart from the beginning.\n\n"
       main  
    fi
}

function change_pythia_version()
{
    local line_numbers line_number current_pythia_version_no_dots new_pythia_version_no_dots
    current_pythia_version_no_dots=${current_pythia_version/./}
    new_pythia_version_no_dots=${new_pythia_version/./}

    for file_examined in $(find ${smash_sources_root_dir} -depth -type f -not \
        \( -path '*/.git/*' -o -name 'CHANGELOG.md' -o -name '*.h' -o -name '*.cc' -o -path "*/3rdparty/*" -o -path "*/bin/*" \));
    do
        line_numbers=$(grep -n -E "[pP]ythia.*${current_pythia_version}" ${file_examined} | cut -d':' -f1);
        if [ -n "${line_numbers}" ]; then
            for line_number in ${line_numbers}
            do
                sed -i "${line_number}s/${current_pythia_version}/${new_pythia_version}/g" ${file_examined};
            done
        fi;
        line_numbers=$(grep -n -E "[pP]ythia.*${current_pythia_version_no_dots}" ${file_examined} | cut -d':' -f1);
        if [ -n "${line_numbers}" ]; then
            for line_number in ${line_numbers}
            do
                sed -i "${line_number}s/${current_pythia_version_no_dots}/${new_pythia_version_no_dots}/g" ${file_examined};
            done
        fi;
    done
}

#=============================================================================
# Auxiliary functions (i.e. not directly called by the main function)

function pretty_exit()
{
    printf "\n \e[1;91mEXITING THE SCRIPT\e[22m $*\e[0m\n" 1>&2
    exit 1
}

function ask_yes_or_no()
{
    local answer go_on
    go_on="y"
    while [ ${go_on} == "y" ]
    do
        printf " \e[1;32my (YES, default option)\e[0m OR \e[1;31mn (NO) \e[0m? (\e[1;36me \e[0mfor \e[1;36mexit\e[0m)\n"
        read answer
        answer=${answer:-y}
        if ! [[ "${answer}" =~ ^[yne]$ ]]
        then
            printf "\e[1;95mInvalid choice!\e[0m\n\n"
        else
            go_on="n"
        fi
    done
    if [ "${answer}" == "e" ]
    then
        pretty_exit
    elif [ "${answer}" == "y" ]
    then
        return 0
    else
        return 1
    fi    
}

function ask_user_smash_source_tree_path()
{
    local tmp_smash_path
    printf "Please, insert the \e[4mabsolute path\e[0m of your SMASH sources. (\e[1;36me \e[0mfor \e[1;36mexit\e[0m)\n"
    read tmp_smash_path
    if [ "${tmp_smash_path}" == "e" ]
    then
        pretty_exit
    fi
    if ! [ -d "${tmp_smash_path}" ]
    then
        printf "Sorry, but the directory that you have entered does not exist. Please, check and try again.\n\n"
        ask_user_smash_source_tree_path
    else
        if ! [ -f "${tmp_smash_path}/src/smash.cc" ]
        then
            printf "Sorry, but probably you did not write the correct path.\n"
            printf "The source directory should contain the file 'src/smash.cc', but the directory you provided does not!\n"
            printf "Please, check and try again.\n\n"
            ask_user_smash_source_tree_path
        else
            smash_sources_root_dir=${tmp_smash_path}
        fi
    fi
}

function reminder_update_changelog()
{
    printf "Please, remind to annotate the update to Pythia %s in the changelog file %s.\n\n" \
    "${new_pythia_version}" "${smash_sources_root_dir}/CHANGELOG.md"
}

function print_closing_message()
{
    printf "\n\e[35mPythia version changed. \e[1;35mPlease, remember to mention this operation in the file CHANGELOG.md!\e[0m\n\n"
}
#===================================================================
# Main function
main
