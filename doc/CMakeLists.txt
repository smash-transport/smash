find_program(DOXYGEN doxygen)
if(DOXYGEN)
   execute_process(COMMAND ${DOXYGEN} --version
      OUTPUT_VARIABLE DOXYGEN_VERSION
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
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
   configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
      ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
      @ONLY
      )
   configure_file(${CMAKE_CURRENT_SOURCE_DIR}/DoxyUser.in
      ${CMAKE_CURRENT_BINARY_DIR}/DoxyUser
      @ONLY
      )
   add_custom_target(doc ${DOXYGEN} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      COMMENT "Building Doxygen Code Documentation"
      VERBATIM
      )
   find_program(SED sed)
   execute_process(COMMAND ${SED} --version OUTPUT_VARIABLE _sed_version)
   if(NOT "${_sed_version}" MATCHES "GNU sed")
      unset(SED CACHE)
      set(SED "SED-NOTFOUND")
      message(WARNING "GNU sed is required to build the user guide")
   endif()
   find_program(GREP grep)
   if(SED)
      file(GLOB files_with_input "${PROJECT_SOURCE_DIR}/src/include/*.h" "${PROJECT_SOURCE_DIR}/src/*.cc"
                                 "${PROJECT_SOURCE_DIR}/AUTHORS" "${PROJECT_SOURCE_DIR}/INSTALL")
      add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/userguide.dox
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/generateuserdocs.sh ${SED}
         ${CMAKE_CURRENT_BINARY_DIR}/userguide.dox ${files_with_input}
         DEPENDS ${files_with_input}
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
         COMMENT "Generating User Guide input file"
         VERBATIM
         )
      add_custom_target(user ${DOXYGEN} ${CMAKE_CURRENT_BINARY_DIR}/DoxyUser
         DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/userguide.dox
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
         COMMENT "Building User Guide"
         VERBATIM
         )
   endif()

   find_program(XARGS xargs)
   find_program(PRINTF printf)
   find_program(SORT sort)
   if(GREP AND XARGS AND PRINTF AND SORT)
      set(DOXYGEN_EXTRACT_ALL "NO")
      set(HAVE_DOT "NO")
      configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
         ${CMAKE_CURRENT_BINARY_DIR}/DoxyfileOnlyDocumented
         @ONLY
         )

      add_custom_target(undocumented ${DOXYGEN}
         ${CMAKE_CURRENT_BINARY_DIR}/DoxyfileOnlyDocumented 2>&1 |
         ${SORT} -u |
         ${GREP} -i warning: | ${GREP} -iv "einhard"
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
         COMMENT "Building Doxygen Code Documentation (only with documented entities) and output doxgen warnings"
         VERBATIM
         )

      add_custom_target(undocumented_count ${DOXYGEN}
         ${CMAKE_CURRENT_BINARY_DIR}/DoxyfileOnlyDocumented 2>&1 |
         ${GREP} -v 3rdparty |
         ${SORT} -u |
         ${GREP} -i warning: | ${GREP} -icv "einhard" |
         ${XARGS} ${PRINTF} "============ There are %d doxygen warnings about missing documentation ============\\n"
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
         COMMENT "Building Doxygen Code Documentation (only with documented entities) and count doxgen warnings"
         VERBATIM
         )
   endif()
   mark_as_advanced(XARGS PRINTF SORT SED GREP DOT_EXECUTABLE DOXYGEN)
else()
   message(STATUS "Doxygen not found.")
endif()
