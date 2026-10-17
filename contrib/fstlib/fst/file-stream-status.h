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
#ifndef FST_FILE_STREAM_STATUS_H_
#define FST_FILE_STREAM_STATUS_H_

#include <cerrno>
#include <cstring>
#include <fstream>
#include <iostream>

#include <fst/status/status.h>
#include <fstream>

namespace fst {

inline ::fst::Status GetFileStreamStatus(const std::ios& stream) {
  if (stream.fail()) {
    return ::fst::InternalError(std::strerror(errno));
  }
  return ::fst::OkStatus();
}

inline ::fst::Status CloseFileStream(std::ifstream& stream) {
  stream.close();
  if (stream.fail()) {
    return ::fst::InternalError(std::strerror(errno));
  }
  return ::fst::OkStatus();
}

inline ::fst::Status CloseFileStream(std::ofstream& stream) {
  stream.close();
  if (stream.fail()) {
    return ::fst::InternalError(std::strerror(errno));
  }
  return ::fst::OkStatus();
}

inline ::fst::Status CloseFileStream(std::fstream& stream) {
  stream.close();
  if (stream.fail()) {
    return ::fst::InternalError(std::strerror(errno));
  }
  return ::fst::OkStatus();
}

}  // namespace fst

#endif  // FST_FILE_STREAM_STATUS_H_
