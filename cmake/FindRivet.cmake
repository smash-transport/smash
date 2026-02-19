########################################################
#
#    Copyright (c) 2022,2025
#      SMASH Team
#
#    BSD 3-clause license
#
########################################################

# cmake-format: off
#=======================================================
# - Try to find Rivet
#
# This module defines:
#
#  RIVET_CONFIG          Path to rivet-config executable
#  RIVET_VERSION         Rivet version string
#  RIVET_INCLUDE_DIR     Rivet include directory
#  RIVET_LIBRARIES       Linker flags from rivet-config --libs
#  Rivet_FOUND           Whether Rivet was found
#
# It also defines the imported target Rivet::Rivet which provides:
#   - SYSTEM include directories
#   - INTERFACE link libraries and flags from rivet-config
#=======================================================
# cmake-format: on

find_program(RIVET_CONFIG rivet-config DOC "Path to the Rivet configuration script")

if(RIVET_CONFIG)
    execute_process(COMMAND ${RIVET_CONFIG} --version OUTPUT_VARIABLE RIVET_VERSION
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${RIVET_CONFIG} --includedir OUTPUT_VARIABLE RIVET_INCLUDE_DIR
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${RIVET_CONFIG} --libs OUTPUT_VARIABLE RIVET_LIBRARIES
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Rivet REQUIRED_VARS RIVET_CONFIG RIVET_INCLUDE_DIR RIVET_LIBRARIES
                                  VERSION_VAR RIVET_VERSION)

mark_as_advanced(RIVET_CONFIG)

# ~~~
# We want now to create a modern Rivet::Rivet target. Since we get a list of linker flags from
# rivet-config, we need to separate out the library names from other flags. Therefore, we
# process the RIVET_LIBRARIES list in the following way:
#  - collect all -L paths to search for libraries;
#  - resolve all -l flags to absolute library paths (using the collected -L paths), and
#  - resolve the -lRivet flag separately to an absolute path to the main Rivet library.
# ~~~
if(Rivet_FOUND AND NOT TARGET Rivet::Rivet)
    # Convert whitespace-separated flags into a proper CMake list. This ensures consistency with how
    # CMake find-modules usually return lists of libraries and include directories.
    separate_arguments(RIVET_LIBRARIES NATIVE_COMMAND "${RIVET_LIBRARIES}")
    # Prepare lists to hold absolute library paths and other linker options
    set(RIVET_LIBRARIES_ABS)
    set(RIVET_LINK_OPTIONS)
    set(RIVET_LIB_SEARCH_DIRS)
    set(RIVET_MAIN_LIB "NOTFOUND") # Will hold -lRivet resolved path, NOTFOUND value is necessary,
                                   # otherwise find_library will not perform the search
    foreach(flag IN LISTS RIVET_LIBRARIES)
        if(flag MATCHES "^-L(.+)")
            list(APPEND RIVET_LIB_SEARCH_DIRS "${CMAKE_MATCH_1}")
            list(REMOVE_DUPLICATES RIVET_LIB_SEARCH_DIRS)
        elseif(flag MATCHES "^-l(.+)")
            set(LIB_NAME "${CMAKE_MATCH_1}")
            if(LIB_NAME STREQUAL "Rivet")
                find_library(RIVET_MAIN_LIB
                             NAMES Rivet
                             PATHS ${RIVET_LIB_SEARCH_DIRS}
                             NO_DEFAULT_PATH)
            else()
                # NO_CACHE is crucial to ensure that find_library searches the specified paths and
                # doesn't return a cached result from a previous search
                unset(LIB_FOUND) # Clear any previous value of LIB_FOUND before calling find_library
                find_library(LIB_FOUND
                             NAMES "${LIB_NAME}"
                             PATHS ${RIVET_LIB_SEARCH_DIRS}
                             NO_DEFAULT_PATH NO_CACHE)
                if(LIB_FOUND)
                    list(APPEND RIVET_LIBRARIES_ABS "${LIB_FOUND}")
                else()
                    message(DEBUG "Rivet library flag -l${LIB_NAME} not found "
                                  "in paths: ${RIVET_LIB_SEARCH_DIRS}"
                                  " -> handing it over to linker as link flag.")
                    # Fallback: pass the -l flag as raw option
                    list(APPEND RIVET_LINK_OPTIONS "${flag}")
                endif()
            endif()
        else()
            # Anything else (e.g., -Wl,-rpath,/path or -framework ...)
            list(APPEND RIVET_LINK_OPTIONS "${flag}")
        endif()
    endforeach()
    if(NOT RIVET_MAIN_LIB)
        message(FATAL_ERROR "Cannot create Rivet::Rivet target: main library not found!")
    endif()
    # Create imported target
    add_library(Rivet::Rivet UNKNOWN IMPORTED)
    set_target_properties(Rivet::Rivet PROPERTIES IMPORTED_LOCATION "${RIVET_MAIN_LIB}")
    if(RIVET_LIBRARIES_ABS)
        set_target_properties(Rivet::Rivet PROPERTIES INTERFACE_LINK_LIBRARIES
                                                      "${RIVET_LIBRARIES_ABS}")
    endif()
    if(RIVET_LINK_OPTIONS)
        set_target_properties(Rivet::Rivet PROPERTIES INTERFACE_LINK_OPTIONS
                                                      "${RIVET_LINK_OPTIONS}")
    endif()
    target_include_directories(Rivet::Rivet SYSTEM INTERFACE ${RIVET_INCLUDE_DIR})
endif()
