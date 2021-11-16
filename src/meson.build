# Copyright (c) 2020-2021 Reed A. Cartwright <reed@cartwright.ht>
# Copyright (c) 2021 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# Process version information
semver = meson.project_version().split('.')
version_src = configure_file(input : 'version.h.in.in',
  output : 'version.h.in',
  configuration : {
    'version_major' : semver[0],
    'version_minor' : semver[1],
    'version_patch' : semver[2],
    'version_int'   : (semver[0].to_int()*1000+ semver[1].to_int())*10000+ \
                       semver[2].to_int(),
    'version_vcs'   : '@VCS_TAG@'

  })

version_file = vcs_tag(input : version_src, output : 'version.h')

inc = include_directories(['.', 'include'])

eigen_dep = dependency('eigen3',
	fallback : ['eigen', 'eigen_dep'],
	version : '>=3.3',
	required : true
)

boost_dep = dependency('boost',
	version : '>=1.47.0',
	required : true
)

json_dep = dependency('json',
  fallback : ['json', 'nlohmann_json_dep'],
  version : '>=3.2.0',
  required : true
)

subdir('lib')

progs=['alignpair', 'genseed', 'msa', 'version', 'sample', 'format']

foreach p : progs
  exe = executable('coati-@0@'.format(p), ['coati-@0@.cc'.format(p), version_file],
    link_with : [libcoati],
    include_directories : inc,
    dependencies : [libcoati_deps],
    cpp_args : ['-DDOCTEST_CONFIG_DISABLE'],
    install : true,
    install_dir : get_option('libexecdir')
  )
endforeach

coati_src = configure_file(input : 'coati.cc.in',
  output : 'coati.cpp',
  configuration : {
    'source_dir' : meson.current_source_dir(),
    'prefix'     : get_option('prefix'),
    'libexecdir' : get_option('libexecdir')
  }
)

executable('coati', [coati_src, version_file],
  include_directories : inc,
  link_with : [libcoati],
  install : true
)
