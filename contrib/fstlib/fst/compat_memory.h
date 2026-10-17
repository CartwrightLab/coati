// Copyright 2005-2024 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the 'License');
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an 'AS IS' BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
#ifndef FST_COMPAT_MEMORY_H_
#define FST_COMPAT_MEMORY_H_

#include <stddef.h>

#include <memory>
#include <type_traits>

namespace fst {

// Defines make_unique_for_overwrite using a standard definition that should be
// compatible with the C++20 definition. That is, all compiling uses of
// `std::make_unique_for_overwrite` should have the same result with
// `fst::make_unique_for_overwrite`. Note that the reverse doesn't
// necessarily hold.
// TODO(kbg): Remove these once we migrate to C++20.

template <typename T>
std::unique_ptr<T> make_unique_for_overwrite() {
  return std::unique_ptr<T>(new T);
}

template <typename T>
std::unique_ptr<T> make_unique_for_overwrite(size_t n) {
  return std::unique_ptr<T>(new std::remove_extent_t<T>[n]);
}

template <typename T>
std::unique_ptr<T> WrapUnique(T* ptr) {
  return std::unique_ptr<T>(ptr);
}

}  // namespace fst

#endif  // FST_COMPAT_MEMORY_H_
