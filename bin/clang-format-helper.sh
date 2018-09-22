#!/bin/bash

fail() {
  echo "ERROR: $*" >&2
  exit 1
}

type clang-format >/dev/null 2>&1 || fail "Clang-format not found"

VER_NUM=`echo $(clang-format -version) | sed -ne "s/[^0-9]*\(\([0-9]\.\)\{0,4\}[0-9][^.]\).*/\1/p"`
if [ "$VER_NUM" != "6.0.0 " ]; then
    fail "Wrong clang-format version: $VER_NUM (6.0.0 is required)"
fi;

FILES_TO_FORMAT="../src/*.cc ../src/include/smash/*.h ../src/tests/*.cc ../src/tests/*.h"

if [ "$1" = "-p" ]; then
    echo "Running clang-format ... [perform from within the /bin directory]"
    for i in $FILES_TO_FORMAT;
    do
        clang-format --verbose -i $i;
    done;

elif [ "$1" = "-t" ]; then
    echo "Testing that clang-format does not change the source code in the working directory ..."
    if diff <(cat $FILES_TO_FORMAT) <(clang-format $FILES_TO_FORMAT) > /dev/null; then
        echo "PASS: No changes to source code by clang-format."
    else
        fail "Clang-format was not properly run on latest commit."
    fi

else
cat <<EOF
Little helper script that wraps clang-format
Usage: ./clang-format-helper
  -p perform automatic formatting (for developers)
  -t test that automatic formatting has been performed (for CI builds)
Note: Run from within the /bin directory (uses relative paths).
EOF
fi
