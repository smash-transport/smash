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
   find_program(GREP grep)
   file(GLOB files_with_input "${PROJECT_SOURCE_DIR}/src/include/*.h")
   add_custom_target(inputoptions ${SED}
      "\\_^\\s*\\(\\*\\|//\\)\\s*!!USER:Input_, \\_^\\s*\\(\\*\\|//\\)\\s*!!/USER:Input_!d" ${files_with_input} |
      ${GREP} -v USER:Input > ${CMAKE_CURRENT_BINARY_DIR}/inputoptions.dox
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      COMMENT "Preparing input options file"
      VERBATIM
      )
   add_custom_target(outputexamples ${SED}
      "\\_^\\s*\\(\\*\\|//\\)\\s*!!USER:Output_, \\_^ \\* !!/USER:Output_!d" ${PROJECT_SOURCE_DIR}/src/outputroutines.cc |
      ${GREP} -v USER:Input > ${CMAKE_CURRENT_BINARY_DIR}/outputexamples.dox
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      COMMENT "Preparing output examples file"
      VERBATIM
      )
   add_custom_target(user ${DOXYGEN} ${CMAKE_CURRENT_BINARY_DIR}/DoxyUser
      DEPENDS inputoptions outputexamples
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      COMMENT "Building User Guide"
      VERBATIM
      )

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
