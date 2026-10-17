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
// A shim for ::fst::StatusOr.

#ifndef FST_STATUS_STATUSOR_H_
#define FST_STATUS_STATUSOR_H_

#ifdef OPENFST_HAS_ABSL

namespace fst {
using ::fst::StatusOr;
}  // namespace fst
#else
#include <optional>
#include <type_traits>
#include <utility>

#include <fst/status/status.h>

namespace fst {

// A minimal shim for ::fst::StatusOr.
//
// API differences:
// 1. Does not do any safety checks, assumes the API is used correctly, e.g.:
//    - If `value()` is called on a `Status` instance, the result is undefined.
//    - Constructing a StatusOr with an OK status is undefined behavior.
// 2. Some methods, such as `value_or`, `IgnoreError`, etc, are not implemented.
//    Feel free to add them as needed.
template <typename T>
class StatusOr {
 public:
  template <typename U = T,
            std::enable_if_t<std::is_constructible<T, U&&>::value, int> = 0>
  // NOLINTNEXTLINE(google-explicit-constructor)
  StatusOr(U&& u) : value_(std::forward<U>(u)) {}
  template <typename U = T,
            std::enable_if_t<std::is_constructible<T, U&&>::value, int> = 0>
  // NOLINTNEXTLINE(google-explicit-constructor)
  StatusOr(StatusOr<U>&& u) {
    if (u.ok()) {
      value_.emplace(std::forward<StatusOr<U> >(u).value());
    } else {
      status_ = std::forward<StatusOr<U> >(u).status();
    }
  }

  // NOLINTNEXTLINE(google-explicit-constructor)
  StatusOr(const T& value) : value_(value) {}
  // NOLINTNEXTLINE(google-explicit-constructor)
  StatusOr(T&& value) : value_(std::move(value)) {}
  // NOLINTNEXTLINE(google-explicit-constructor)
  StatusOr(const Status& status) : status_(status) {}
  // NOLINTNEXTLINE(google-explicit-constructor)
  StatusOr(Status&& status) : status_(std::move(status)) {}

  explicit StatusOr() = default;
  StatusOr(const StatusOr<T>&) = default;
  StatusOr& operator=(const StatusOr<T>&) = default;
  StatusOr(StatusOr<T>&&) noexcept = default;
  StatusOr& operator=(StatusOr<T>&&) noexcept = default;

  bool ok() const { return value_.has_value(); }

  const T& value() const& { return value_.value(); }
  T& value() & { return value_.value(); }
  T value() && { return std::move(value_).value(); }

  const Status& status() const& { return status_; }
  Status status() && { return std::move(status_); }

  const T* operator->() const { return &*value_; }
  T* operator->() { return &*value_; }

  const T& operator*() const& { return *value_; }
  T& operator*() & { return *value_; }
  T operator*() && { return *std::move(value_); }

 private:
  Status status_ = Status(StatusCode::kUnknown, "");
  std::optional<T> value_;
};
}  // namespace fst

#endif

#endif  // FST_STATUS_STATUSOR_H_
