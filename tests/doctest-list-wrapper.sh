#!/usr/bin/env sh

doctest_cmd=$1
doctest_bin=$2
doctest_txt=$3

regex="^[=\[]"

if [ "x${doctest_cmd}" = "xupdate" ]; then
  # Extract tests names and print them to stdout
  "$doctest_bin" --list-test-cases | grep -v "$regex" > "$doctest_txt"
  exit 0
fi

if [ "x$doctest_cmd" = "xtest" ]; then
  # compare existing names to extracted names
  "$doctest_bin" --list-test-cases | grep -v "$regex" | diff -q "${doctest_txt}" - >/dev/null
  res=$?
  if [ $res -ne 0 ]; then
    echo "FAILURE: Unit tests are out of date."
    echo "         Run 'meson compile update-tests' to fix."
    exit $res
  fi
  exit 0
fi

echo "ERROR: Unknown command '${doctest_cmd}'."
exit 1
