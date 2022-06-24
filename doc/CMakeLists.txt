########################################################
#
#    Copyright (c) 2014,2017-2018,2020-2022
#      SMASH Team
#
#    BSD 3-clause license
#
#########################################################

find_program(DOXYGEN doxygen)
if(DOXYGEN)
    execute_process(COMMAND ${DOXYGEN} --version OUTPUT_VARIABLE DOXYGEN_VERSION
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    if("${DOXYGEN_VERSION}" VERSION_LESS "1.8.4")
        message(WARNING "Doxygen ${DOXYGEN_VERSION} found, but at least 1.8.4 is expected.")
    else()
        message(STATUS "Doxygen ${DOXYGEN_VERSION} found at ${DOXYGEN}")
    endif()

    find_program(DOT_EXECUTABLE dot)
    if(DOT_EXECUTABLE)
        set(HAVE_DOT "YES")
    else()
        set(HAVE_DOT "NO")
    endif()

    set(DOXYGEN_EXTRACT_ALL "YES")
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                   @ONLY)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/DoxyUser.in ${CMAKE_CURRENT_BINARY_DIR}/DoxyUser
                   @ONLY)
    add_custom_target(doc
                      ${DOXYGEN} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                      COMMENT "Building Doxygen Code Documentation"
                      VERBATIM)
    find_program(SED sed)
    execute_process(COMMAND ${SED} --version OUTPUT_VARIABLE _sed_version)
    if(NOT "${_sed_version}" MATCHES "GNU sed")
        unset(SED CACHE)
        set(SED "SED-NOTFOUND")
        message(WARNING "GNU sed is required to build the user guide")
    endif()
    find_program(GREP grep)
    if(NOT GREP)
        message(WARNING "grep not found, docs might be incomplete/faulty")
    endif()
    find_program(BASH bash)
    if(NOT BASH)
        message(WARNING "bash not found, docs might be incomplete/faulty")
    endif()
    if(SED)
        include(UserInputFiles.txt)
        add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/userguide.dox
                           COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/generateuserdocs.sh ${SED}
                                   ${CMAKE_CURRENT_BINARY_DIR}/userguide.dox ${files_with_input}
                           DEPENDS ${files_with_input}
                           WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                           COMMENT "Generating User Guide input file"
                           VERBATIM)
        add_custom_target(user
                          ${DOXYGEN} ${CMAKE_CURRENT_BINARY_DIR}/DoxyUser
                          COMMAND ${CMAKE_COMMAND} -E env source_dir=${PROJECT_SOURCE_DIR}
                                  binary_dir=${CMAKE_CURRENT_BINARY_DIR} bash
                                  ${CMAKE_CURRENT_SOURCE_DIR}/set_smash_version.sh
                          DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/userguide.dox
                          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                          COMMENT "Building User Guide"
                          VERBATIM)
    endif()

    find_program(XARGS xargs)
    find_program(PRINTF printf)
    find_program(SORT sort)
    if(GREP
       AND XARGS
       AND PRINTF
       AND SORT)
        set(DOXYGEN_EXTRACT_ALL "NO")
        set(HAVE_DOT "NO")
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
                       ${CMAKE_CURRENT_BINARY_DIR}/DoxyfileOnlyDocumented @ONLY)

        add_custom_target(undocumented
                          ${DOXYGEN}
                          ${CMAKE_CURRENT_BINARY_DIR}/DoxyfileOnlyDocumented
                          2>&1
                          |
                          ${SORT}
                          -u
                          |
                          ${GREP}
                          -i
                          warning:
                          |
                          ${GREP}
                          -iv
                          "einhard"
                          ||
                          true
                          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                          COMMENT "Building Doxygen Code Documentation (only with documented entities) and output doxgen warnings"
                          VERBATIM)

        add_custom_target(undocumented_count
                          ${DOXYGEN}
                          ${CMAKE_CURRENT_BINARY_DIR}/DoxyfileOnlyDocumented
                          2>&1
                          |
                          ${GREP}
                          -v
                          3rdparty
                          |
                          ${SORT}
                          -u
                          |
                          ${GREP}
                          -i
                          warning:
                          |
                          ${GREP}
                          -icv
                          "einhard"
                          |
                          ${XARGS}
                          ${PRINTF}
                          "============ There are %d doxygen warnings about missing documentation ============\\n"
                          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                          COMMENT "Building Doxygen Code Documentation (only with documented entities) and count doxgen warnings"
                          VERBATIM)

        add_custom_target(undocumented_test
                          ${DOXYGEN}
                          ${CMAKE_CURRENT_BINARY_DIR}/DoxyfileOnlyDocumented
                          2>&1
                          |
                          ${GREP}
                          -v
                          3rdparty
                          |
                          ${SORT}
                          -u
                          |
                          ${GREP}
                          -i
                          warning:
                          |
                          ${GREP}
                          -icv
                          "einhard"
                          |
                          ${GREP}
                          -q
                          "0"
                          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                          COMMENT "Testing Doxygen Code Documentation for CI. Errors if there are warnings."
                          VERBATIM)
    endif()
    mark_as_advanced(XARGS
                     PRINTF
                     SORT
                     SED
                     GREP
                     DOT_EXECUTABLE
                     DOXYGEN)
else()
    message(STATUS "Doxygen not found.")
endif()
