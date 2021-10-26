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

gbenchmark_dep = dependency('google-benchmark',
  fallback : ['google-benchmark', 'google_benchmark_dep'],
	version : '>=1.5.2',
	required : get_option('dev_benchmark'))

bm_src = configure_file(
  input : 'benchmark_main.cc.in',
  output : 'benchmark_main.cc',
  # depfile : 'benchmark_main.cc.in',
  configuration : {
    'bm_data' : meson.current_source_dir()+'/data',
  }
)

benchmark_exe = executable('all_benchmarks',
  meson.project_build_root()+'/benchmark/benchmark_main.cc',
  include_directories : inc,
  dependencies : [gbenchmark_dep, libcoati_deps],
  link_with : [libcoati]
  build_by_default : false
)

benchmark_list = files('libcoati-benchmark-tests.txt')

benchmark_list_wrapper = find_program('benchmark-list-wrapper.sh')

run_target('update-benchmark',
  command : [benchmark_list_wrapper, 'update', benchmark_exe, benchmark_list]
)

benchmark('[pretest] Checking benchmark list completeness',
  benchmark_list_wrapper,
  args : ['test', benchmark_exe, benchmark_list],
  priority : 1000
)

fs = import('fs')
benchmarks = fs.read(benchmark_list).strip().split('\n')
foreach benchmark_case : benchmarks
  benchmark(f'[libcoati] @benchmark_case@', benchmark_exe,
    args : [
      f'--benchmark_filter=@benchmark_case@',
      f'--benchmark_out=@benchmark_case@.json'.replace('/','_'),
      '--benchmark_out_format=json'
    ],
    timeout : 600
  )
endforeach
