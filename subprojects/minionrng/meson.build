project('minionrng', 'cpp')

minion_inc = include_directories('.')
minionrng = static_library('minionrng',
  'minion.cc')

minionrng_dep = declare_dependency(include_directories : minion_inc,
  link_with : minionrng)
