#!/bin/sh
SED="$1"
output="$2"
shift 2
$SED -e '/\/\*[*!] *\\Userguide\>/,/\*\//!d' "$@" > $output
