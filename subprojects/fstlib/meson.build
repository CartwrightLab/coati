
project('fstlib', 'cpp',
  version : '1.8.1',
  default_options: ['cpp_std=c++17', 'warning_level=0'])

fst_inc = include_directories('.')
fstlib = static_library('fstlib',
  'compat.cc', 'encode.cc', 'flags.cc', 'fst-types.cc', 'fst.cc',
  'mapped-file.cc', 'properties.cc', 'symbol-table-ops.cc', 'symbol-table.cc',
  'util.cc', 'weight.cc')

fstlib_dep = declare_dependency(include_directories : fst_inc,
  link_with : fstlib)
