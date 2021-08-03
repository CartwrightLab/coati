// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// See www.openfst.org for extensive documentation on this weighted
// finite-state transducer library.
//
// Google-style logging declarations and inline definitions.

#ifndef FST_LIB_LOG_H_
#define FST_LIB_LOG_H_

#include <cassert>
#include <iostream>
#include <string>

#include <fst/types.h>
#include <fst/flags.h>

DECLARE_int32(v);

class LogMessage {
 public:
  LogMessage(const std::string &type) : fatal_(type == "FATAL") {
    std::cerr << type << ": ";
  }
  ~LogMessage() {
    std::cerr << std::endl;
    if(fatal_)
      exit(1);
  }
  std::ostream &stream() { return std::cerr; }

 private:
  bool fatal_;
};

#define LOG(type) LogMessage(#type).stream()
#define VLOG(level) if ((level) <= FST_FLAGS_v) LOG(INFO)

// Checks
inline void FstCheck(bool x, const char* expr,
                const char *file, int line) {
  if (!x) {
    LOG(FATAL) << "Check failed: \"" << expr
               << "\" file: " << file
               << " line: " << line;
  }
}

#define FSTCHECK(x) FstCheck(static_cast<bool>(x), #x, __FILE__, __LINE__)
#define FSTCHECK_EQ(x, y) FSTCHECK((x) == (y))
#define FSTCHECK_LT(x, y) FSTCHECK((x) < (y))
#define FSTCHECK_GT(x, y) FSTCHECK((x) > (y))
#define FSTCHECK_LE(x, y) FSTCHECK((x) <= (y))
#define FSTCHECK_GE(x, y) FSTCHECK((x) >= (y))
#define FSTCHECK_NE(x, y) FSTCHECK((x) != (y))

// Debug checks
#define FSTDCHECK(x) assert(x)
#define FSTDCHECK_EQ(x, y) FSTDCHECK((x) == (y))
#define FSTDCHECK_LT(x, y) FSTDCHECK((x) < (y))
#define FSTDCHECK_GT(x, y) FSTDCHECK((x) > (y))
#define FSTDCHECK_LE(x, y) FSTDCHECK((x) <= (y))
#define FSTDCHECK_GE(x, y) FSTDCHECK((x) >= (y))
#define FSTDCHECK_NE(x, y) FSTDCHECK((x) != (y))


// Ports
#define ATTRIBUTE_DEPRECATED __attribute__((deprecated))

#endif  // FST_LIB_LOG_H_
