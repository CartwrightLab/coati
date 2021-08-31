#!/bin/#

# find all tests and add them to a single file

dir=$1
keyword=$2

# all_files="$(ls "${dir}"/*.cc)"
all_files=$(find "${dir}"/ -name "*.cc" -exec basename {} \;)

# find all tests
for cpp_file in ${all_files[*]}
do
    file=$(basename "${cpp_file}")
    cat "${dir}/${file}" | awk ' BEGIN { state = 0; last = ""; }
        $0 ~ /^'"${keyword}"'\(/ { print last; state = 1; }
          { if (state == 1) print; }
        $0 ~ /^}/ { if (state) state = 2; }
          { last = $0; }
        '
done
