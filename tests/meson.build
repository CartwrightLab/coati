# Copyright (c) 2021 Reed A. Cartwright <reed@cartwright.ht>
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

doctest_exe = executable('libcoati-doctest', ['libcoati-doctest.cc', version_file, libcoati_sources],
  include_directories : inc,
  dependencies : [doctest_dep, libcoati_deps],
  build_by_default : false
)

test_list = files('libcoati-doctest-tests.txt')

doctest_list_wrapper = find_program('doctest-list-wrapper.sh')

run_target('update-tests',
  command : [doctest_list_wrapper, 'update', doctest_exe, test_list]
)

test('[pretest] Checking test list completeness', doctest_list_wrapper,
  args : ['test', doctest_exe, test_list],
  is_parallel : false,
  priority : 1000
)

fs = import('fs')
tests = fs.read(test_list).strip().split('\n')
foreach test_case : tests
  test(f'[libcoati] @test_case@',
    doctest_exe,
    args : [f'--test-case=@test_case@'],
    is_parallel : false
  )
endforeach
