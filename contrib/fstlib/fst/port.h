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
// This header file defines a set of macros for checking the presence of
// compiler and platform features. Such macros are used to produce portable
// code based on the presence or lack of a given feature.

#ifndef FST_PORT_H_
#define FST_PORT_H_

#if defined(__clang__)
#define FST_LTO_VISIBILITY_PUBLIC [[clang::lto_visibility_public]]
#else
#define FST_LTO_VISIBILITY_PUBLIC
#endif

#endif  // FST_PORT_H_
