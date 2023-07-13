/*
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
*/

#ifndef SEMIRING_HPP
#define SEMIRING_HPP

#include "utils.hpp"

namespace coati::semiring {
// Define semiring types and their operations.

/**
 * @brief Linear semiring class.
 *
 * |         Set          | + | * | + neutral element | * neutral element |
 * | Reals U {-INF, +INF} | + | * |         0         |         1         |
 *
 */
class linear {
   public:
    static constexpr float plus(float x, float y) { return x + y; }
    static constexpr float plus(float x, float y, float z) { return x + y + z; }
    static constexpr float times(float x) { return x; }
    template <typename... Args>
    static constexpr float times(float first, Args... args) {
        return first * times(args...);
    }
    static constexpr float power(float x, float y) { return std::pow(x, y); }
    static constexpr float power(float x, size_t y) {
        return std::pow(x, static_cast<float>(y));
    }
    static constexpr float zero() { return 0.0f; }
    static constexpr float one() { return 1.0f; }
};
/**
 * @brief Log semiring class.
 *
 * |         Set          |     +     | * | + neutral element | * neutral elem |
 * | Reals U {-INF, +INF} | logSumExp | + |        -INF       |         0      |
 *
 */
class log {
   public:
    static float plus(float x, float y) {
        return coati::utils::log_sum_exp(x, y);
    }
    static float plus(float x, float y, float z) {
        return coati::utils::log_sum_exp(coati::utils::log_sum_exp(x, y), z);
    }
    static constexpr float times(float x) { return x; }
    template <typename... Args>
    static constexpr float times(float first, Args... args) {
        return first + times(args...);
    }
    static constexpr float power(float x, float y) { return x * y; }
    static constexpr float power(float x, size_t y) {
        return x * static_cast<float>(y);
    }
    static constexpr float zero() {
        return std::numeric_limits<float_t>::lowest();
    }
    static constexpr float one() { return 0.0f; }

    static float from_linearf(float x) { return ::logf(x); }
    static float to_linearf(float x) { return ::expf(x); }
    static float from_linear_1mf(float x) { return std::log1pf(-x); }
    static float to_linear_1mf(float x) { return -std::expm1f(x); }
};
/**
 * @brief Tropical semiring class.
 *
 * |         Set          |  +  | * | + neutral element | * neutral element |
 * | Reals U {-INF, +INF} | max | + |        -INF       |         0      |
 *
 */
class tropical {
   public:
    static constexpr float plus(float x, float y) { return std::max(x, y); }
    static constexpr float plus(float x, float y, float z) {
        return std::max(std::max(x, y), z);
    }
    static constexpr float times(float x) { return x; }
    template <typename... Args>
    static constexpr float times(float first, Args... args) {
        return first + times(args...);
    }
    static constexpr float power(float x, float y) { return x * y; }
    static constexpr float power(float x, size_t y) {
        return x * static_cast<float>(y);
    }
    static constexpr float zero() {
        return std::numeric_limits<float_t>::lowest();
    }
    static constexpr float one() { return 0.0f; }

    static float from_linearf(float x) { return ::logf(x); }
    static float to_linearf(float x) { return ::expf(x); }
    static float from_linear_1mf(float x) { return std::log1pf(-x); }
    static float to_linear_1mf(float x) { return -std::expm1f(x); }
};
}  // namespace coati::semiring

#endif
