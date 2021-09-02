# Copyright (c) 2020 Reed A. Cartwright <reed@cartwright.ht>
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

inc = include_directories('include')

# version configure_file can't be run on a different dir (conf), see:
# https://github.com/mesonbuild/meson/issues/30
configure_file(input : 'version.h.in',
  output : 'version.h',
  configuration : version_conf)

conf_inc = include_directories('.')

subdir('lib')

coati_version = executable('coati-version',
  'coati-version.cc',
  include_directories : [inc, conf_inc],
  link_with : libcoati)

coati_alignpair = executable('coati-alignpair',
  'coati-alignpair.cc',
  include_directories : [inc, conf_inc],
  link_with : libcoati,
  dependencies : [doctest_dep, eigen_dep, boost_dep, fstlib_dep],
  install : true)

coati_msa = executable('coati-msa',
  'coati-msa.cc',
  include_directories : [inc, conf_inc],
  link_with : libcoati,
  dependencies : [doctest_dep, eigen_dep, boost_dep, fstlib_dep],
  install : true)

# coati main script
coati_conf = configuration_data()
coati_conf.set('MESON_CURRENT_SOURCE_DIR', meson.current_source_dir())
coati_conf.set('MESON_INSTALL_PREFIX', get_option('prefix'))
coati_conf.set('MESON_INSTALL_LIBEXECDIR', get_option('libexecdir'))

configure_file(input : meson.current_source_dir() + '/coati.cc.in',
  output : 'coati.cc',
  configuration : coati_conf)

coati = executable('coati',
  meson.current_build_dir() + '/coati.cc',
  include_directories : [inc, conf_inc],
  link_with : libcoati,
  dependencies : [doctest_dep, eigen_dep, boost_dep],
  install : true)

configure_file(input : meson.current_source_dir() + '/coati-doctest.cc.in',
  output : 'coati-doctest.cc',
  configuration : tests_conf,
  depfile : 'libcoati_dep')

coati_doctest = executable('coati-doctest',
  meson.project_build_root() + '/src/coati-doctest.cc',
  include_directories : inc,
  link_with : libcoati,
  dependencies : [doctest_dep, eigen_dep, boost_dep, fstlib_dep])

test('libcoati', coati_doctest)
