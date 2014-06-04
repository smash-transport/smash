# - Try to find ROOT
# This will define:
#
#  ROOT_FOUND
#  ROOT_CONFIG
#  ROOT_INCLUDE_DIR
#  ROOT_LIBRARIES
#  ROOT_VERSION

#=============================================================================
# Copyright © 2014  SMASH Team
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
find_program(ROOT_CONFIG NAMES root-config HINTS $ENV{ROOTSYS}/bin)
if(ROOT_CONFIG)
   execute_process(COMMAND "${ROOT_CONFIG}" --prefix
      OUTPUT_VARIABLE ROOTSYS OUTPUT_STRIP_TRAILING_WHITESPACE)

   execute_process(COMMAND "${ROOT_CONFIG}" --version
      OUTPUT_VARIABLE ROOT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
   string(REPLACE "/" "." ROOT_VERSION "${ROOT_VERSION}")

   execute_process(COMMAND "${ROOT_CONFIG}" --incdir
      OUTPUT_VARIABLE ROOT_INCDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
   set(ROOT_INCLUDE_DIR ${ROOT_INCDIR} CACHE PATH "path to ROOT includes")

   execute_process(COMMAND "${ROOT_CONFIG}" --glibs
      OUTPUT_VARIABLE ROOT_LINK_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
   set(ROOT_LIBRARIES ${ROOT_LINK_FLAGS} CACHE PATH "linker flags for ROOT")
endif()

FIND_PACKAGE_HANDLE_STANDARD_ARGS(ROOT
   REQUIRED_VARS ROOT_CONFIG ROOT_INCLUDE_DIR ROOT_LIBRARIES
   VERSION_VAR ROOT_VERSION
   FAIL_MESSAGE "required ROOT library not found"
   )

mark_as_advanced(ROOT_CONFIG)
