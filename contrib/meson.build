doctest_dep = declare_dependency(include_directories: include_directories('doctest', is_system : true))

librandom = static_library(
    'random', ['random/random.cpp'],
    include_directories : include_directories('random')
)

librandom_dep = declare_dependency(
    include_directories : include_directories('random'),
    link_with : librandom
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

cli_dep = declare_dependency(include_directories : include_directories('CLI', is_system : true))
