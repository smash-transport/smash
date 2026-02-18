########################################################
#
#    Copyright (c) 2026
#      SMASH Team
#
#    BSD 3-clause license
#
#########################################################

include("${CMAKE_CURRENT_LIST_DIR}/CheckCCompilerFlag.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/CheckCXXCompilerFlag.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/CheckMicCCompilerFlag.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/CheckMicCXXCompilerFlag.cmake")

# Tests whether a compiler flag is supported for the given language and compiler family. Does not
# modify any targets or global compiler flags. The result is returned in the variable specified by
# RESULT. Syntax:
# ~~~
#    check_compiler_flag_is_supported(LANGUAGE <C|CXX>
#                                     FLAG     <flag>
#                                     RESULT   <out-var>
#                                     [COMPILER <AUTO|MIC>])
# ~~~
function(check_compiler_flag_is_supported)
    cmake_parse_arguments(arg_of
                          ""
                          "LANGUAGE;FLAG;RESULT;COMPILER"
                          ""
                          ${ARGN})

    if(NOT arg_of_LANGUAGE)
        message(FATAL_ERROR "check_compiler_flag_is_supported: LANGUAGE is required (C or CXX)")
    else()
        string(TOUPPER "${arg_of_LANGUAGE}" _lang)
        if(NOT _lang STREQUAL "C" AND NOT _lang STREQUAL "CXX")
            message(FATAL_ERROR "check_compiler_flag_is_supported: LANGUAGE must be C or CXX")
        endif()
        if(NOT CMAKE_${_lang}_COMPILER_LOADED)
            message(FATAL_ERROR "check_compiler_flag_is_supported: ${_lang} language is not enabled"
            )
        endif()
    endif()

    if(NOT arg_of_FLAG)
        message(FATAL_ERROR "check_compiler_flag_is_supported: FLAG is required")
    endif()

    if(NOT arg_of_RESULT)
        message(FATAL_ERROR "check_compiler_flag_is_supported: RESULT variable name is required")
    endif()

    if(arg_of_COMPILER)
        string(TOUPPER "${arg_of_COMPILER}" _compiler)
    else()
        set(_compiler "AUTO")
    endif()

    if(NOT _compiler STREQUAL "AUTO" AND NOT _compiler STREQUAL "MIC")
        message(FATAL_ERROR "check_compiler_flag_is_supported: COMPILER must be AUTO or MIC")
    endif()

    string(MAKE_C_IDENTIFIER "SMASH_HAVE_${_lang}_${arg_of_FLAG}" _cache_var)

    if(NOT DEFINED CACHE{${_cache_var}})
        if(_compiler STREQUAL "MIC")
            if(_lang STREQUAL "C")
                check_mic_c_compiler_flag("${arg_of_FLAG}" ${_cache_var})
            else()
                check_mic_cxx_compiler_flag("${arg_of_FLAG}" ${_cache_var})
            endif()
        else()
            if(_lang STREQUAL "C")
                check_c_compiler_flag("${arg_of_FLAG}" ${_cache_var})
            else()
                check_cxx_compiler_flag("${arg_of_FLAG}" ${_cache_var})
            endif()
        endif()
        # Although this might be done by try_compile inside the called macros above, it is harmless
        # to enforce it here and in this way it is obvious that this function is caching the result.
        set(${_cache_var} ${${_cache_var}} CACHE BOOL
                                                 "Whether compiler supports flag ${arg_of_FLAG}")
        if(NOT ${_cache_var})
            message(ATTENTION "Flag ${arg_of_FLAG} not supported by ${_lang} compiler."
                              " It will not be used.")
        endif()
    endif()

    if(${_cache_var})
        set(${arg_of_RESULT} TRUE PARENT_SCOPE)
    else()
        set(${arg_of_RESULT} FALSE PARENT_SCOPE)
    endif()
endfunction()
