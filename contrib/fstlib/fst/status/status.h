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
// See www.openfst.org for extensive documentation on this weighted
// finite-state transducer library.
//
// A shim for ::fst::Status.

#ifndef FST_STATUS_STATUS_H_
#define FST_STATUS_STATUS_H_

#ifdef OPENFST_HAS_ABSL

namespace fst {
using ::absl::AbortedError;
using ::absl::AlreadyExistsError;
using ::absl::CancelledError;
using ::absl::DataLossError;
using ::absl::DeadlineExceededError;
using ::absl::FailedPreconditionError;
using ::absl::InternalError;
using ::absl::InvalidArgumentError;
using ::absl::NotFoundError;
using ::absl::OkStatus;
using ::absl::OutOfRangeError;
using ::absl::PermissionDeniedError;
using ::absl::ResourceExhaustedError;
using ::absl::Status;
using ::absl::StatusCode;
using ::absl::UnauthenticatedError;
using ::absl::UnavailableError;
using ::absl::UnimplementedError;
using ::absl::UnknownError;
}  // namespace fst

#else

#include <iostream>
#include <string>

#include <string_view>

namespace fst {

enum class StatusCode : int {
  kOk = 0,
  kCancelled = 1,
  kUnknown = 2,
  kInvalidArgument = 3,
  kDeadlineExceeded = 4,
  kNotFound = 5,
  kAlreadyExists = 6,
  kPermissionDenied = 7,
  kResourceExhausted = 8,
  kFailedPrecondition = 9,
  kAborted = 10,
  kOutOfRange = 11,
  kUnimplemented = 12,
  kInternal = 13,
  kUnavailable = 14,
  kDataLoss = 15,
  kUnauthenticated = 16,
  kDoNotUseReservedForFutureExpansionUseDefaultInSwitchInstead_ = 20
};

class Status final {
 public:
  Status() = default;
  Status(const Status&) = default;
  Status(Status&&) noexcept = default;
  Status& operator=(const Status&) = default;
  Status& operator=(Status&&) noexcept = default;

  Status(StatusCode code, std::string_view msg)
      : code_(code), message_(std::string(msg)) {}

  bool ok() const { return code_ == StatusCode::kOk; }
  StatusCode code() const { return code_; }
  std::string_view message() const { return message_; }
  std::string ToString() const { return std::string(message_); }

  friend bool operator==(const Status& lhs, const Status& rhs) {
    return lhs.code_ == rhs.code_ && lhs.message_ == rhs.message_;
  }
  friend bool operator!=(const Status& lhs, const Status& rhs) {
    return !(lhs == rhs);
  }

 private:
  StatusCode code_ = StatusCode::kOk;
  std::string message_;
};

inline std::ostream& operator<<(std::ostream& os, const Status& status) {
  return os << status.ToString();
}

inline Status OkStatus() { return Status(); }
inline Status AbortedError(std::string_view message) {
  return Status(StatusCode::kAborted, message);
}
inline Status AlreadyExistsError(std::string_view message) {
  return Status(StatusCode::kAlreadyExists, message);
}
inline Status CancelledError(std::string_view message) {
  return Status(StatusCode::kCancelled, message);
}
inline Status DataLossError(std::string_view message) {
  return Status(StatusCode::kDataLoss, message);
}
inline Status DeadlineExceededError(std::string_view message) {
  return Status(StatusCode::kDeadlineExceeded, message);
}
inline Status FailedPreconditionError(std::string_view message) {
  return Status(StatusCode::kFailedPrecondition, message);
}
inline Status InternalError(std::string_view message) {
  return Status(StatusCode::kInternal, message);
}
inline Status InvalidArgumentError(std::string_view message) {
  return Status(StatusCode::kInvalidArgument, message);
}
inline Status NotFoundError(std::string_view message) {
  return Status(StatusCode::kNotFound, message);
}
inline Status OutOfRangeError(std::string_view message) {
  return Status(StatusCode::kOutOfRange, message);
}
inline Status PermissionDeniedError(std::string_view message) {
  return Status(StatusCode::kPermissionDenied, message);
}
inline Status ResourceExhaustedError(std::string_view message) {
  return Status(StatusCode::kResourceExhausted, message);
}
inline Status UnavailableError(std::string_view message) {
  return Status(StatusCode::kUnavailable, message);
}
inline Status UnimplementedError(std::string_view message) {
  return Status(StatusCode::kUnimplemented, message);
}
inline Status UnknownError(std::string_view message) {
  return Status(StatusCode::kUnknown, message);
}

}  // namespace fst
#endif  // OPENFST_HAS_ABSL

#endif  // FST_STATUS_STATUS_H_
