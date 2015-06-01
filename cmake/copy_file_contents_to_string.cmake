########################################################
#
#    Copyright (c) 2014-2015
#      SMASH Team
#
#    BSD 3-clause license
# 
######################################################### 

file(READ "${INPUT_FILE}" TEXT)
string(REPLACE ";" "\\;" TEXT "${TEXT}")
string(REPLACE "\n" ";" TEXT "${TEXT}")
set(output "// This file was generated from ${INPUT_FILE}. Do not modify.\n")
set(output "${output}const char ${NAME}[] =")
foreach(line ${TEXT})
   set(output "${output}\n  \"${line}\\n\"")
endforeach()
file(WRITE "${OUTPUT_FILE}" "${output};\n")
