# OpenFst's core library, libfst, built from the release tarball
# subprojects/openfst.wrap downloads: the sources of OpenFst's own
# src/lib/Makefile.am, as a static library.
project('openfst', 'cpp',
  version : '1.8.5',
  license : 'Apache-2.0',
  meson_version : '>=0.58.0',
  default_options : ['cpp_std=c++17'])

fst_lib = static_library('fst',
  ['src/lib/compat.cc', 'src/lib/encode.cc', 'src/lib/flags.cc',
   'src/lib/fst-types.cc', 'src/lib/fst.cc', 'src/lib/mapped-file.cc',
   'src/lib/properties.cc', 'src/lib/symbol-table-ops.cc',
   'src/lib/symbol-table.cc', 'src/lib/util.cc', 'src/lib/weight.cc'],
  include_directories : include_directories('src/include'),
  override_options : ['warning_level=0']
)

fst_dep = declare_dependency(
  include_directories : include_directories('src/include', is_system : true),
  link_with : fst_lib
)
