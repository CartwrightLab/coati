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

#  COdon-aware Alignment Transducer Project Configuration
# Init Project
project('COATI', 'cpp',
	version : '0.1',
	license : 'MIT',
	meson_version : '>=0.56.0',
	default_options: [ 'buildtype=debugoptimized', 'cpp_std=c++17' ])

#
#
# pretest = executable('pretest', configure?,
#   build_by_default : true,
# 	build_always_stale : true)

# Find Boost installation with library program_options
# https://mesonbuild.com/Dependencies.html#boost
boost_dep = dependency('boost',
	version : '>=1.47.0',
	modules : ['program_options'],
	required : true)

# Find Eigen library
eigen_dep = dependency('eigen3',
	fallback : ['eigen', 'eigen_dep'],
	version : '>=3.3',
	required : true)

doctest_dep = dependency('doctest',
	fallback : ['doctest', 'doctest_dep'],
	version : '>=2.4.6',
	required : true)

fstlib_proj = subproject('fstlib')
fstlib_dep = fstlib_proj.get_variable('fstlib_dep')
fstlib_dep = fstlib_dep.as_system('system') # TODO: this didn't fix the warnings

minionrng_proj = subproject('minionrng')
minionrng_dep = minionrng_proj.get_variable('minionrng_dep')

subdir('conf')
subdir('src')

#subdir('tests')
#subdir('doc')
