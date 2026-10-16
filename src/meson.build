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
    'version_str'   : '@VCS_TAG@'
  })

version_file = vcs_tag(input : version_src, output : 'version.h')

inc = include_directories(['.', 'include'])

eigen_dep = dependency('eigen3',
	fallback : ['eigen', 'eigen_dep'],
	version : '>=3.3',
	required : true,
  include_type : 'system'
)

boost_dep = dependency('boost',
	version : '>=1.67.0',
	required : true,
  include_type : 'system'
)

json_dep = dependency('json',
  fallback : ['json', 'nlohmann_json_dep'],
  version : '>=3.2.0',
  required : true,
  include_type : 'system'  
)

# OpenFst: an installed one when it works, else OpenFst 1.8.5 downloaded and
# built by Meson (subprojects/openfst.wrap). OpenFst installs no pkg-config
# or CMake file, so an installed library is used only if a program using
# OpenFst the way COATi does links against it. 1.8.2 and 1.8.3 do not: their
# fst/fst.h assigns a raw pointer to a std::unique_ptr, which current
# compilers reject.
openfst_check = '''
#include <fst/fstlib.h>

int main() {
  fst::StdVectorFst a;
  const auto s = a.AddState();
  const auto t = a.AddState();
  a.SetStart(s);
  a.AddArc(s, fst::StdArc(1, 1, 0.5, t));
  a.SetFinal(t, 0.0);
  fst::EncodeMapper<fst::StdArc> encoder(fst::kEncodeLabels);
  fst::Encode(&a, &encoder);
  fst::Decode(&a, encoder);
  fst::StdVectorFst b;
  fst::Compose(a, a, &b);
  fst::StdVectorFst path;
  fst::ShortestPath(b, &path);
  fst::SymbolTable symbols;
  symbols.AddSymbol("A");
  return fst::Verify(path) && symbols.Find("A") == 0 ? 0 : 1;
}
'''
cxx = meson.get_compiler('cpp')
openfst_opt = get_option('system_openfst')
fst_dep = dependency('', required : false)
if not openfst_opt.disabled()
  fst_lib = cxx.find_library('fst', has_headers : ['fst/fstlib.h'],
    required : openfst_opt)
  if fst_lib.found() and cxx.links(openfst_check, dependencies : fst_lib,
      name : 'installed OpenFst works with COATi')
    fst_dep = fst_lib
  elif openfst_opt.enabled()
    error('system_openfst=enabled, but no installed OpenFst (1.8.4 or later) links a test program')
  endif
endif
if not fst_dep.found()
  fst_dep = subproject('openfst').get_variable('fst_dep')
endif

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
