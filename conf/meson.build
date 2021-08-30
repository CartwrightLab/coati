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

################################################################################
# version

version_list = meson.project_version().split('.')
major = version_list[0].to_int()
minor = version_list[1].to_int()

version_conf = configuration_data()
version_conf.set('COATI_SEMVER_MAJOR', major)
version_conf.set('COATI_SEMVER_MINOR', minor)

if version_list.length() < 3
  # prerelease
  version_patch = run_command('git', 'describe', '--dirty', '--always')
  version_int = (major * 1000 + minor) * 10000
  version_str = meson.project_version() + \
    '.0-prerelease+g' + \
    version_patch.stdout().strip()

  version_conf.set('COATI_SEMVER_PATCH', 0)
  version_conf.set('COATI_SEMVER_INTEGER', version_int)
  version_conf.set('COATI_SEMVER_VERSION', version_str)
else
  patch = version_list[2].to_int()
  version_int = ((major * 1000 + minor) * 1000 + patch) * 10
  version_conf.set('COATI_SEMVER_PATCH', patch)
  version_conf.set('COATI_SEMVER_INTEGER', version_int)
  version_conf.set('COATI_SEMVER_VERSION', meson.project_version())
endif

# `configure_file` can't be run here yet, see:
#  https://github.com/mesonbuild/meson/issues/30
