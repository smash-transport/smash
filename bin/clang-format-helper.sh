#!/bin/sh

fail() {
  echo "$*" >&2
  exit 1
}

type clang-format >/dev/null 2>&1 || fail "Clang-format not found"

VER_NUM=`echo $(clang-format -version) | sed -ne "s/[^0-9]*\(\([0-9]\.\)\{0,4\}[0-9][^.]\).*/\1/p"`
if [ "$VER_NUM" != "6.0.0 " ]; then
    fail "Wrong clang-format version: $VER_NUM (6.0.0 is required)" 
fi;

run_cf() {
    for i in ../src/*.cc ../src/include/smash/*.h ../src/tests/*.cc ../src/tests/*.h;
    do 
        clang-format -verbose -i $i; 
    done;
}

if [ "$1" == "" ] || [ "$1" == "-h" ]; then
cat <<EOF
Little helper script that wraps clang-format
usage: ./clang-format-helper 
  -p perform automatic formatting (for developers)
  -t test that automatic formatting has been performed (for CI)
  -h display this help 
EOF
elif [ "$1" == "-p" ]; then
    echo "Running clang-format ... [perform only from wihtin the /bin directory]"
    run_cf
elif [ "$1" == "-t" ]; then
    echo "Testing that clang-format does no changes to source code ..."
    run_cf
    if git diff --exit-code --name-status; then
        echo "No changes to source code by clang-format. PASS."
    else
        fail "Clang-format was not properly run on latest commit."
    fi
fi
