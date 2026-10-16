doctest_dep = declare_dependency(include_directories: include_directories('doctest', is_system : true))

librandom = static_library(
    'random', ['random/random.cpp'],
    include_directories : include_directories('random')
)

librandom_dep = declare_dependency(
    include_directories : include_directories('random'),
    link_with : librandom
)

cli_dep = declare_dependency(include_directories : include_directories('CLI', is_system : true))
