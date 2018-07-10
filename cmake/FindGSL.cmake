# - Try to find GSL (GNU Scientific Library)
# This will define:
#
#  GSL_FOUND
#  GSL_INCLUDE_DIR
#  GSL_LIBRARY
#  GSL_CBLAS_LIBRARY
#
#  This module is based on the module provided by cmake starting with version 3.2
#  See https://cmake.org/licensing for copyright info.


#=============================================================================
# Copyright © 2018  SMASH Team
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#  * Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
#  * Neither the names of contributing organizations nor the
#    names of its contributors may be used to endorse or promote products
#    derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS ``AS IS''
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#=============================================================================


include(FindPackageHandleStandardArgs)

# first check if GSL_ROOT_DIR is set (either as environment variable or
# supplied as cmake-option. If so, use it
if (EXISTS "$ENV{GSL_ROOT_DIR}")
  file( TO_CMAKE_PATH "$ENV{GSL_ROOT_DIR}" GSL_ROOT_DIR )
  set( GSL_ROOT_DIR "${GSL_ROOT_DIR}" CACHE PATH "Prefix for GSL installation")

 elseif (EXISTS "${GSL_ROOT_DIR}" )
  file( TO_CMAKE_PATH ${GSL_ROOT_DIR} GSL_ROOT_DIR )
  set( GSL_ROOT_DIR "${GSL_ROOT_DIR}" CACHE PATH "Prefix for GSL installation")
endif()

# no user supplied path. Try to find gsl on our own.
if ( NOT EXISTS "${GSL_ROOT_DIR}" )
  set( GSL_USE_PKGCONFIG ON )
endif()


if (GSL_USE_PKGCONFIG)
  find_package(PkgConfig)
  pkg_check_modules( GSL gsl )
  if (EXISTS "${GSL_INCLUDE_DIR}")
    get_filename_component( GSL_ROOT_DIR "${GSL_INCLUDE_DIR}" PATH CACHE )
  endif()
endif()

find_path( GSL_INCLUDE_DIR
  NAMES gsl
  HINTS ${GSL_ROOT_DIR}/include ${GSL_INCLUDE_DIR}
  )


find_library( GSL_LIBRARY
  NAMES gsld gsl
  HINTS ${GSL_ROOT_DIR}/lib ${GSL_LIBDIR}
  )

find_library( GSL_CBLAS_LIBRARY
  NAMES gslcblas cblas
  HINTS ${GSL_ROOT_DIR}/lib ${GSL_LIBDIR}
  )

set (GSL_INCLUDE_DIRS ${GSL_INCLUDE_DIR} )
set (GSL_LIBRARIES ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} )

# we could probably also search in the filepath for the version.
if (NOT GSL_VERSION)
  find_program( GSL_CONFIG_EXE
    NAMES gsl-config
    HINTS ${GSL_ROOT_DIR}/bin ${GSL_ROOT_DIR}
    )
  if (EXISTS ${GSL_CONFIG_EXE})
    execute_process(
      COMMAND "${GSL_CONFIG_EXE}" "--version"
      RESULT_VARIABLE GSL_PROC_STATUS
      OUTPUT_VARIABLE GSL_VERSION
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
   endif()
endif()

FIND_PACKAGE_HANDLE_STANDARD_ARGS(GSL
  FOUND_VAR GSL_FOUND
  REQUIRED_VARS
    GSL_INCLUDE_DIR
    GSL_LIBRARY
    GSL_CBLAS_LIBRARY
  VERSION_VAR
    GSL_VERSION
   )

 mark_as_advanced(GSL_CBLAS_LIBRARY GSL_CONFIG_EXE GSL_INCLUDE_DIRS GSL_LIBRARIES )

