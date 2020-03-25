include_directories(${CMAKE_SOURCE_DIR})
include(CheckCXXCompilerFlag)

set(all "ALL")
if(DEFINED ENV{APPVEYOR})
   set(all)
endif()

function(vir_add_run_target target)
   if("${CMAKE_GENERATOR}" MATCHES "Visual Studio")
      # do nothing. This just clutters the solution explorer
   else()
      add_custom_target(run_${target} ${all}
         ${CMAKE_CTEST_COMMAND} -V -R ${target}
         #${target} -v
         DEPENDS ${target}
         COMMENT "Execute ${target} test"
         VERBATIM
         )
   endif()
endfunction()

function(vir_apply_flags name std)
   if(NOT MSVC)
      if(CMAKE_CXX_COMPILER MATCHES "/(icpc|icc)$")
         set_target_properties(${name} PROPERTIES COMPILE_FLAGS "-std=${std} -Wall -Wextra -fp-model source")
      else()
         set_target_properties(${name} PROPERTIES COMPILE_FLAGS "-std=${std} -Wall -Wextra")
      endif()
   endif()
endfunction()

function(vir_add_test name)
  add_executable(${name} ${name}.cpp)
  vir_apply_flags(${name} "c++11")
  add_test(NAME ${name} COMMAND ${name} -v)
  vir_add_run_target(${name})

  if(NOT MSVC)
    check_cxx_compiler_flag("-std=c++14" supports14)
    if(supports14)
      add_executable(${name}-14 ${name}.cpp)
      vir_apply_flags(${name}-14 "c++14")
      add_test(NAME ${name}-14 COMMAND ${name}-14 -v)
      vir_add_run_target(${name}-14)
    endif()
  endif()
endfunction()

add_executable(plotdist plotdist.cpp)
vir_apply_flags(plotdist "c++11")
add_test(NAME plotdist
   COMMAND ${CMAKE_COMMAND}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/plotdist.cmake
   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
vir_add_run_target(plotdist)

vir_add_test(checks)
vir_add_test(empty)
