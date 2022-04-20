########################################################
#
#    Copyright (c) 2022
#      SMASH Team
#
#    BSD 3-clause license
# 
#########################################################

# Redefine CMake message to have the possibility to suppress
# informational messages (not warnings, errors or custom attentions)
function(message)
    if(NOT WIN32)
        string(ASCII 27 Esc)
        set(DEFAULT_COL "${Esc}[m")
        set(YELLOW      "${Esc}[93m")
    endif()
    list(GET ARGV 0 _message_type)
    if(_message_type MATCHES "^(ATTENTION|(FATAL|SEND)_ERROR|(AUTHOR_)?WARNING|DEPRECATION|NOTICE|STATUS|VERBOSE|DEBUG|TRACE|CHECK_(START|PASS|FAIL))$")
        list(REMOVE_AT ARGV 0)
        set(_message_start_at 1)
    else()
        set(_message_start_at 0)
        unset(_message_type)
    endif()
    # Do as CMake does, from documentation of message command:
    #   "If more than one message string is given,
    #   they are concatenated into a single message
    #   with no separator between the strings."
    # To stay general is not trivial, since one could
    # pass strings containing semicolons here and we do not want
    # to "loose" them treating as separators in lists
    # https://cmake.org/cmake/help/latest/manual/cmake-language.7.html#lists
    math(EXPR _stop_at "${ARGC}-1")
    foreach(index RANGE ${_message_start_at} ${_stop_at})
        set(_text_of_message "${_text_of_message}${ARGV${index}}")
    endforeach()
    if (NOT MESSAGE_QUIET OR "${_message_type}" MATCHES "(ERROR|WARNING|ATTENTION)")
        if("${_message_type}" STREQUAL "ATTENTION")
            _message(STATUS "${YELLOW}${_text_of_message}${DEFAULT_COL}")
        else()
            _message(${_message_type} "${_text_of_message}")
        endif()
    endif()
endfunction()

# Add utility function to add a compiler flag in a sound way
# (i.e. testing if it is supported) and possibly warn or fail if
# it is not. Syntax:
#
#    add_compiler_flags_if_supported(<flag(s)> [VERBOSE] [ON_FAILURE <value>]
#                                    [C_FLAGS <var>] [CXX_FLAGS <var>])
#
# Passing either C_FLAGS or CXX_FLAGS or both make the function
# add the flag to the passed flag variable(s). If none is passed
# the flag is added to both CMAKE_C_FLAGS and CMAKE_CXX_FLAGS.
# ON_FAILURE accepted values are QUIET|WARN|FATAL and WARN is the
# used one if nothing is passed.
#
# TECHNICAL NOTES:
#  1. The function is prepared to work for flags containing ';' as well,
#     but this case is excluded at the moment, since the AddCompilerFlag
#     macro, which is used internally, does not support such possibility.
#     It was decided to postpone a fix till a real need occurs.
#     The PARSE_ARGV variant of cmake_parse_arguments can only be used
#     in functions and not in macro, hence this must be a function.
#  2. This function has variables semantically very similar to the
#     AddCompilerFlag macro. However, the same name cannot/should not
#     be used, because there would be a clash due to how macros work.
#     Functions, instead, have their own scope, but changing the
#     AddCompilerFlag macro into a function would also not naively
#     work, since its <...>_RESULT variables have to be set in the
#     calling scope. This is possible using set(... PARENT_SCOPE)
#     but was not done. As convention, local variables here have been
#     prefixed with a double underscore "__".
get_filename_component(_currentDir "${CMAKE_CURRENT_LIST_FILE}" PATH)
include("${_currentDir}/AddCompilerFlag.cmake")
function(add_compiler_flags_if_supported)
    # Parse arguments and do logic to set up needed variables
    cmake_parse_arguments(PARSE_ARGV 0 "_" "VERBOSE" "ON_FAILURE;C_FLAGS;CXX_FLAGS" "")
    list(LENGTH __UNPARSED_ARGUMENTS __number_of_flags)
    if(__number_of_flags EQUAL 0)
        message(FATAL_ERROR "No flag passed to add_compiler_flags_if_supported!")
    else()
        set(__flags "${__UNPARSED_ARGUMENTS}")
    endif()
    if(__C_FLAGS)
        set(__c_flags "${__C_FLAGS}")
    endif()
    if(__CXX_FLAGS)
        set(__cxx_flags "${__CXX_FLAGS}")
    endif()
    if(NOT __C_FLAGS AND NOT __CXX_FLAGS)
        set(__c_flags "CMAKE_C_FLAGS")
        set(__cxx_flags "CMAKE_CXX_FLAGS")
    endif()
    if(__ON_FAILURE)
        set(__mode "${__ON_FAILURE}")
    else()
        set(__mode "WARN")
    endif()
    if("${__mode}" STREQUAL "QUIET")
        unset(__action_on_failure)
    elseif("${__mode}" STREQUAL "WARN")
        set(__action_on_failure "ATTENTION")
        set(__unused_flag_message " and this will not be used")
    elseif("${__mode}" STREQUAL "FATAL")
        set(__action_on_failure "FATAL_ERROR")
    else()
        message(FATAL_ERROR "Syntax error for add_compiler_flags_if_supported (wrong verbosity)")
    endif()
    # Finally check/add flags and report to user
    foreach(__flag ${__flags})
        if(__flag MATCHES ";")
            if(DEFINED __action_on_failure)
                message(ATTENTION "Compiler flags containing a semicolon cannot be added, ignoring '${__flag}'.")
            endif()
            continue()
        endif()
        set(MESSAGE_QUIET ON)
        if(DEFINED __c_flags AND DEFINED __cxx_flags)
            AddCompilerFlag("${__flag}" C_FLAGS ${__c_flags} CXX_FLAGS ${__cxx_flags}
                                       C_RESULT __c_result  CXX_RESULT __cxx_result)
        elseif(DEFINED __c_flags)
            AddCompilerFlag("${__flag}" C_FLAGS ${__c_flags} C_RESULT __c_result)
        elseif(DEFINED __cxx_flags)
            AddCompilerFlag("${__flag}" CXX_FLAGS ${__cxx_flags} CXX_RESULT __cxx_result)
        else()
            message(FATAL_ERROR "Unexpected case for add_compiler_flags_if_supported")
        endif()
        unset(MESSAGE_QUIET)
        if(DEFINED __c_flags)
            if(__c_result)
                set(__added_to "${__c_flags}")
            else()
                set(__unsupported_lang "C")
            endif()
        endif()
        if(DEFINED __cxx_flags)
            if(__cxx_result)
                set(__added_to ${__added_to} "${__cxx_flags}")
            else()
                set(__unsupported_lang ${__unsupported_lang} "C++")
            endif()
        endif()
        if(DEFINED __action_on_failure)
            if(DEFINED __unsupported_lang)
                string(REPLACE ";" "/" __unsupported_lang "${__unsupported_lang}") # list(JOIN ...) available from CMake 3.12
                message(${__action_on_failure} "Your ${__unsupported_lang} compiler does not support the '${__flag}' flag${__unused_flag_message}!")
            endif()
        endif()
        if(__VERBOSE)
            if(DEFINED __added_to)
                string(REPLACE ";" " and " __added_to "${__added_to}") # list(JOIN ...) available from CMake 3.12
                message(STATUS "Compiler flag '${__flag}' added to ${__added_to}.")
            endif()
        endif()
        # Unset variables needed for next iteration
        unset(__c_result)
        unset(__cxx_result)
        unset(__added_to)
        unset(__unsupported_lang)
    endforeach()
    if(DEFINED __c_flags)
        set(${__c_flags} "${${__c_flags}}" PARENT_SCOPE)
    endif()
    if(DEFINED __cxx_flags)
        set(${__cxx_flags} "${${__cxx_flags}}" PARENT_SCOPE)
    endif()
endfunction()
