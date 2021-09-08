doctest_dep = declare_dependency(include_directories: include_directories('doctest'))

minionrng = static_library(
    'minionrng', ['minionrng/minion.cc'],
    include_directories : include_directories('minionrng')
)

minionrng_dep = declare_dependency(
    include_directories : include_directories('minionrng'),
    link_with : minionrng
)

fstlib = static_library(
    'fstlib', ['fstlib/compat.cc', 'fstlib/encode.cc', 'fstlib/flags.cc',
    'fstlib/fst-types.cc', 'fstlib/fst.cc', 'fstlib/mapped-file.cc',
    'fstlib/properties.cc', 'fstlib/symbol-table-ops.cc',
    'fstlib/symbol-table.cc', 'fstlib/util.cc', 'fstlib/weight.cc'],
  include_directories : include_directories('fstlib'),
  override_options : 'warning_level=0'
)

fstlib_dep = declare_dependency(
  include_directories : include_directories('fstlib'),
  link_with : fstlib
)

fstlib_dep = fstlib_dep.as_system('system')
