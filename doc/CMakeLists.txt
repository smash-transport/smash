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
   configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
      ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
      @ONLY
      )

   add_custom_target(dox ${DOXYGEN} Doxyfile
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Building Doxygen Code Documentation"
      VERBATIM
      )
else()
   message(STATUS "Doxygen not found.")
endif()
