#!/usr/bin/env sh

benchmark_cmd=$1
benchmark_bin=$2
benchmark_txt=$3

if [ "x${benchmark_cmd}" = "xupdate" ]; then
    # Extract tests names and print them to stdout
    "$benchmark_bin" --benchmark_list_tests=true > "$benchmark_txt"
    exit 0
fi

if [ "x$benchmark_cmd" = "xtest" ]; then
    # Compare existing names to extracted names
    "$benchmark_bin" --benchmark_list_tests=true | diff -q "${benchmark_txt}" - >/dev/null
    res=$?
    if [ $res -ne 0 ]; then
        echo "FAILURE: Benchmark tests are out of date."
        echo "         Run 'meson compile update-benchmark' to fix."
        exit $res
    fi
    exit 0
fi

echo "ERROR: Unknown command '${benchmark_cmd}'."
exit 1
