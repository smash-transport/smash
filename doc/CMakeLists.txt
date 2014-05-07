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
   set(DOXYGEN_EXTRACT_ALL "YES")
   configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
      ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
      @ONLY
      )
   add_custom_target(doc ${DOXYGEN} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      COMMENT "Building Doxygen Code Documentation"
      VERBATIM
      )

   find_program(GREP grep)
   find_program(XARGS xargs)
   find_program(PRINTF printf)
   find_program(SORT sort)
   if(GREP AND XARGS AND PRINTF AND SORT)
      set(DOXYGEN_EXTRACT_ALL "NO")
      configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
         ${CMAKE_CURRENT_BINARY_DIR}/DoxyfileOnlyDocumented
         @ONLY
         )

      add_custom_target(undocumented ${DOXYGEN}
         ${CMAKE_CURRENT_BINARY_DIR}/DoxyfileOnlyDocumented 2>&1 |
         ${SORT} -u |
         ${GREP} -ic warning: |
         ${XARGS} ${PRINTF} "============ There are %d undocumented entities ============\\n"
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
         COMMENT "Building Doxygen Code Documentation (only documented entities)"
         VERBATIM
         )
   endif()
else()
   message(STATUS "Doxygen not found.")
endif()
