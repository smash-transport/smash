########################################################
#
#    Copyright (c) 2022,2024-2026
#      SMASH Team
#
#    BSD 3-clause license
#
########################################################

# Redefine CMake message to have the possibility to suppress informational messages (not warnings,
# errors or custom attentions)
function(message)
    if(NOT WIN32)
        string(ASCII 27 Esc)
        set(DEFAULT_COL "${Esc}[m")
        set(YELLOW "${Esc}[93m")
    endif()
    list(GET ARGV 0 _message_type)
    if(_message_type
       MATCHES
       "^(ATTENTION|(FATAL|SEND)_ERROR|(AUTHOR_)?WARNING|DEPRECATION|NOTICE|STATUS|VERBOSE|DEBUG|TRACE|CHECK_(START|PASS|FAIL))$"
    )
        list(REMOVE_AT ARGV 0)
        set(_message_start_at 1)
    elseif(_message_type MATCHES "^CONFIGURE_LOG$")
        set(MESSAGE_QUIET ON)
    else()
        set(_message_start_at 0)
        unset(_message_type)
    endif()
    # Do as CMake does, from documentation of message command: "If more than one message string is
    # given, they are concatenated into a single message with no separator between the strings." To
    # stay general is not trivial, since one could pass strings containing semicolons here and we do
    # not want to "lose" them treating as separators in lists
    # https://cmake.org/cmake/help/latest/manual/cmake-language.7.html#lists
    math(EXPR _stop_at "${ARGC}-1")
    foreach(index RANGE ${_message_start_at} ${_stop_at})
        set(_text_of_message "${_text_of_message}${ARGV${index}}")
    endforeach()
    if(NOT MESSAGE_QUIET OR "${_message_type}" MATCHES "(ERROR|WARNING|ATTENTION)")
        if("${_message_type}" STREQUAL "ATTENTION")
            _message(STATUS "${YELLOW}${_text_of_message}${DEFAULT_COL}")
        else()
            _message(${_message_type} "${_text_of_message}")
        endif()
    endif()
    if(_message_type MATCHES "^CONFIGURE_LOG$")
        set(MESSAGE_QUIET OFF)
    endif()
endfunction()

# Add function to add a compiler flag to a target only if the compiler supports it. Syntax:
# ~~~
# target_add_compiler_flag_if_supported(TARGETS   <target> [<target> ...]   # One or more targets
#                                       LANGUAGE <C|CXX>                    # Language of the flags
#                                       FLAGS    <flag> [<flag> ...]        # Flags to test and add
#                                       [RESULT <variable>]                 # Variable filled if specified
#                                       [COMPILER <AUTO|MIC>]               # Default AUTO
#                                       [SCOPE <PRIVATE|INTERFACE|PUBLIC>]) # Default PRIVATE
# For each target and each flag:
# 1. If the target already has the flag, skip and print a message
# 2. Else, check if the compiler supports the flag using check_compiler_flag_is_supported()
# 3. If supported, add the flag to the target with target_compile_options()
# 4. SCOPE is applied to all targets uniformly
# 5. If RESULT is specified, this variable is set to a list of TRUE/FALSE values indicating whether
#    the flag was added or not to each target. The length of this list is the number of targets passed
#    multiplied by the number of flags passed. The order is the same as the order of targets and flags
#    passed. A flag skipped because already present is considered as added, hence TRUE in the result list.
# ~~~
include("${CMAKE_CURRENT_LIST_DIR}/CheckCompilerFlag.cmake")
function(target_add_compiler_flag_if_supported)
    cmake_parse_arguments(arg_of
                          ""
                          "LANGUAGE;COMPILER;SCOPE;RESULT"
                          "TARGETS;FLAGS"
                          ${ARGN})

    if(NOT arg_of_TARGETS)
        message(FATAL_ERROR "TARGETS is required")
    endif()
    if(NOT arg_of_LANGUAGE)
        message(FATAL_ERROR "LANGUAGE is required (C or CXX)")
    endif()
    if(NOT arg_of_FLAGS)
        message(FATAL_ERROR "FLAGS must contain at least one flag")
    else()
        foreach(_flag IN LISTS arg_of_FLAGS)
            if(_flag MATCHES "\\$<")
                message(FATAL_ERROR " \n"
                                    " Generator expressions are not allowed in FLAGS:\n   ${_flag}\n"
                                    " Pass raw compiler flags only.\n")
            endif()
        endforeach()
        # If the same flag is passed more than once, just ignore duplicates
        list(REMOVE_DUPLICATES arg_of_FLAGS)
    endif()

    if(arg_of_SCOPE)
        string(TOUPPER "${arg_of_SCOPE}" _scope)
    else()
        set(_scope PRIVATE)
    endif()
    if(NOT _scope MATCHES "^(PRIVATE|INTERFACE|PUBLIC)$")
        message(FATAL_ERROR "SCOPE must be PRIVATE, PUBLIC, or INTERFACE")
    endif()

    if(arg_of_COMPILER)
        set(_compiler "${arg_of_COMPILER}")
    else()
        set(_compiler AUTO)
    endif()

    if(arg_of_RESULT)
        set(_result_list "")
    endif()

    foreach(_target IN LISTS arg_of_TARGETS)
        # Collect existing target compile options and interface compile options (if needed) to check
        # if the flag is already present and avoid adding it twice. Note that possibly existing
        # generator expressions in target options remain unevaluated strings here, but this is fine
        # as this function is intended to be used at configure time.
        set(_existing_opts "")
        get_target_property(_options ${_target} COMPILE_OPTIONS)
        if(NOT _options STREQUAL "_options-NOTFOUND")
            list(APPEND _existing_opts ${_options})
        endif()
        if(_scope MATCHES "^(PUBLIC|INTERFACE)$")
            get_target_property(_options ${_target} INTERFACE_COMPILE_OPTIONS)
            if(NOT _options STREQUAL "_options-NOTFOUND")
                list(APPEND _existing_opts ${_options})
            endif()
        endif()

        foreach(_flag IN LISTS arg_of_FLAGS)
            list(FIND _existing_opts "${_flag}" _found_index)
            if(_found_index GREATER -1)
                message(ATTENTION "The target '${_target}' already has the flag '${_flag}',"
                                  " skipping it.")
                if(arg_of_RESULT)
                    list(APPEND _result_list TRUE) # Skipped flags considered “added”
                endif()
                continue()
            endif()
            set(MESSAGE_QUIET ON)
            check_compiler_flag_is_supported(LANGUAGE ${arg_of_LANGUAGE}
                                             FLAG "${_flag}"
                                             RESULT _flag_supported
                                             COMPILER ${_compiler})
            unset(MESSAGE_QUIET)
            if(NOT _flag_supported)
                if(arg_of_RESULT)
                    list(APPEND _result_list FALSE)
                endif()
                continue()
            else()
                # Note that the generator expression is expanded as late as possible and the
                # following one is meant to avoid passing C-only options to C++ compiler when
                # compiling C++ sources and vice versa. For pure C or C++ targets it is irrelevant,
                # but for mixed one it would be wrong not do do so.
                target_compile_options(${_target} ${_scope}
                                       $<$<COMPILE_LANGUAGE:${arg_of_LANGUAGE}>:${_flag}>)
                if(arg_of_RESULT)
                    list(APPEND _result_list TRUE)
                endif()
            endif()
        endforeach()
    endforeach()

    # If RESULT variable is specified, set it in the calling scope
    if(arg_of_RESULT)
        set(${arg_of_RESULT} "${_result_list}" PARENT_SCOPE)
    endif()
endfunction()
