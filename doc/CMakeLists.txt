########################################################
#
#    Copyright (c) 2014,2017-2018,2020-2025
#      SMASH Team
#
#    BSD 3-clause license
#
########################################################

# NOTE: If the minimum required version of Doxygen is increased, the documentation configuration
# files must be updated. To do so, run 'doxygen -u Doxyfile.in' and 'doxygen -u DoxyUser.in' from
# the doc folder (in the top-level directory of the smash repository), making sure to be using the
# lowest required Doxygen version.
find_package(Doxygen 1.9.1 OPTIONAL_COMPONENTS dot)
if(DOXYGEN_FOUND)
    if(DOXYGEN_VERSION VERSION_GREATER "1.9.5")
        message(ATTENTION "Your Doxygen version has not been tested. "
                          "Only versions 1.9.[1-2] are known to work properly. "
                          "Warnings and/or unexpected behavior might occur.")
    elseif(DOXYGEN_VERSION VERSION_GREATER "1.9.2")
        message(ATTENTION "Doxygen version ${DOXYGEN_VERSION} is known "
                          "not to build SMASH documentation fully as expected. "
                          "Only versions 1.9.[1-2] are known to work properly. "
                          "Use at your own risk.")
    endif()

    # Store Doxygen executable full path in a variable for later
    get_target_property(DOXYGEN Doxygen::doxygen IMPORTED_LOCATION)

    # In the following, we set CMake variables that seems to be unused, but which instead are used
    # by configure_file to set-up Doxy* files. E.g. HAVE_DOT, DOT_EXECUTABLE, DOXYGEN_EXTRACT_ALL
    if(Doxygen_dot_FOUND)
        set(HAVE_DOT "YES")
        get_target_property(DOT_EXECUTABLE Doxygen::dot IMPORTED_LOCATION)
    else()
        set(HAVE_DOT "NO")
    endif()
    set(DOXYGEN_EXTRACT_ALL "YES")
    # Prepare list of files to be used as input for user guide and development documentation
    set(COMMON_GENERIC_FILES
        "README.md"
        "LICENSE.md"
        "INSTALL.md"
        "CHANGELOG.md"
        "CONTRIBUTING.md"
        "AUTHORS.md"
        "doc/index.dox")
    set(COMMON_TECHNICAL_FILES "containers/README.md" "bin/benchmarks/README.md")
    list(APPEND
         DOCUMENTATION_FILES_DEVEL
         "doc/mainpage.dox"
         "${COMMON_GENERIC_FILES}"
         "src/include/smash/"
         "3rdparty/einhard/einhard.hpp"
         "src/tests/unittest.h"
         "src/tests/setup.h"
         "src/"
         "3rdparty/README.md"
         "${COMMON_TECHNICAL_FILES}")
    # The user guide has the main page in the doc folder in the userguide.dox file, which has to be
    # configured by CMake in order to possibly set the correct link to the codebase source. Hence
    # here we already include the output of such configuration which will be put in the build/doc
    # folder as userguide_mainpage.dox file. All the /*!\Userguide blocks, instead, are extracted in
    # the build/doc/userguide.dox folder. Both files have to be included, the former first and the
    # second after the md files.
    set(USERGUIDE_MAINPAGE "${CMAKE_CURRENT_BINARY_DIR}/userguide_mainpage.dox")
    list(APPEND
         DOCUMENTATION_FILES_USER
         "${USERGUIDE_MAINPAGE}"
         "${COMMON_GENERIC_FILES}"
         "${CMAKE_CURRENT_BINARY_DIR}/userguide.dox"
         "${COMMON_TECHNICAL_FILES}")
    set(INDENTATION "                         ")
    list(JOIN DOCUMENTATION_FILES_DEVEL "\\\n${INDENTATION}" DOCUMENTATION_FILES_DEVEL)
    list(JOIN DOCUMENTATION_FILES_USER "\\\n${INDENTATION}" DOCUMENTATION_FILES_USER)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                   @ONLY)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/DoxyUser.in ${CMAKE_CURRENT_BINARY_DIR}/DoxyUser
                   @ONLY)
    add_custom_target(doc
                      ${DOXYGEN} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                      COMMENT "Building Doxygen Code Documentation"
                      VERBATIM)

    # Collect requirement tools (sed, bash) and set up user-guide target
    find_program(SED sed)
    execute_process(COMMAND ${SED} --version OUTPUT_VARIABLE _sed_version)
    if(NOT "${_sed_version}" MATCHES "GNU sed")
        unset(SED CACHE)
        set(SED "SED-NOTFOUND")
        message(WARNING "GNU sed is required to build the user guide")
    endif()
    find_program(BASH bash)
    if(NOT BASH)
        message(WARNING "bash not found, docs might be incomplete/faulty")
    endif()
    if(SED)
        # We want that the shipped public user guide contains a link to the repository, but the
        # documentation might also be built in-between stable releases and it might be misleading to
        # leave there the same sentence. Hence fine-tune here.
        if(NOT SMASH_VERSION MATCHES "^SMASH-[0-9]+([.][0-9]+)*$")
            string(CONCAT SENTENCE_WITH_CODEBASE_LINK
                          "\\note This user guide has not been generated "
                          "from a stable version of the codebase.")
            set(SENTENCE_FOR_TARGET_COLOR "--yellow")
            set(SENTENCE_FOR_TARGET "Not on stable release, add note to user guide")
        else()
            set(LINK_TO_GITHUB
                "https://github.com/smash-transport/smash/releases/tag/${SMASH_VERSION}")
            string(CONCAT SENTENCE_WITH_CODEBASE_LINK
                          "The source code for the stable release from "
                          "which this user guide originated can be found "
                          "<a href=\"${LINK_TO_GITHUB}\">here</a>.")
            set(SENTENCE_FOR_TARGET_COLOR "--cyan")
            set(SENTENCE_FOR_TARGET "On stable release, add link to codebase in user guide")
        endif()
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/userguide.dox ${USERGUIDE_MAINPAGE} @ONLY)
        include(UserInputFiles.cmake)
        add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/userguide.dox
                           COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/prepare_userguide_source.bash ${SED}
                                   ${CMAKE_CURRENT_BINARY_DIR}/userguide.dox -- ${files_with_input}
                           DEPENDS ${files_with_input}
                           WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                           COMMENT "Generating User Guide input file"
                           VERBATIM)

        add_custom_target(user
                          ${DOXYGEN} ${CMAKE_CURRENT_BINARY_DIR}/DoxyUser
                          COMMAND ${CMAKE_COMMAND} -E cmake_echo_color
                                  "${SENTENCE_FOR_TARGET_COLOR}" --bold "${SENTENCE_FOR_TARGET}"
                          DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/userguide.dox
                          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                          COMMENT "Building User Guide"
                          VERBATIM)
    endif()

    # Collect requirement tools and set up utility targets
    find_program(GREP grep)
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
    message(ATTENTION "Without Doxygen, it won't be possible to compile the documentation.")
endif()
