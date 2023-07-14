/*
# Copyright (c) 2020-2021 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com>
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

#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cassert>

#include <vector>

namespace coati {

/**
 * @brief Matrix class with custom constructors.
 *
 */
template <class T>
class Matrix {
   public:
    Matrix() = default;
    Matrix(std::size_t rows, std::size_t cols, T value = static_cast<T>(0))
        : rows_(rows), cols_(cols), data_(rows * cols, value) {}

    template< class InputIt >
    Matrix(std::size_t rows, std::size_t cols, InputIt first, InputIt last)
        : rows_(rows), cols_(cols), data_(first, last) {
    }

    // copy constructor
    Matrix(const Matrix&) = default;
    // move constructor
    Matrix(Matrix&&) noexcept = default;
    // assignment operator
    Matrix& operator=(const Matrix&) = default;
    // move assignment operator
    Matrix& operator=(Matrix&&) noexcept = default;

    Matrix(std::initializer_list<std::initializer_list<T>> init_list)
        : rows_(init_list.size()), cols_(init_list.begin()->size()) {
        data_.resize(rows_ * cols_);

        size_t i = 0, j = 0;
        for(const auto& row : init_list) {
            for(const auto& val : row) {
                data_[i * cols_ + j] = val;
                j++;
            }
            i++;
            j = 0;
        }
    }

    // destructor
    ~Matrix() = default;

    T operator()(std::size_t row, std::size_t col) const {
        assert(row < rows_ && col < cols_);
        return data_[row * cols_ + col];
    }
    T& operator()(std::size_t row, std::size_t col) {
        assert(row < rows_ && col < cols_);
        return data_[row * cols_ + col];
    }
    bool operator==(const Matrix& mat) const {
        if((rows_ != mat.rows_) || (cols_ != mat.cols_)) {
            return false;
        }

        for(std::size_t i = 0; i < rows_ * cols_; i++) {
            if(data_[i] != mat.data_[i]) {
                return false;
            }
        }
        return true;
    }

    void resize(std::size_t rows, std::size_t cols,
                T value = static_cast<T>(0)) {
        rows_ = rows;
        cols_ = cols;
        data_.resize(rows * cols);
        data_.assign(data_.size(), value);
    }

    [[nodiscard]] std::size_t rows() const { return rows_; }
    [[nodiscard]] std::size_t cols() const { return cols_; }

   private:
    std::size_t rows_{0}, cols_{0};
    std::vector<T> data_;
};  // class matrix

////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Tensor class - multidimensional matrices.
 *
 */
template <class T>
class Tensor {
   public:
    Tensor(std::size_t dims, std::size_t rows, std::size_t cols, T value = 0.0f)
        : dims_(dims), rows_(rows), cols_(cols) {
        data_.resize(dims_ * rows_ * cols_);
        data_.assign(data_.size(), value);
    }
    // copy constructor
    Tensor(const Tensor& tens)
        : dims_(tens.dims_), rows_(tens.rows_), cols_(tens.cols_) {
        memcpy(&data_[0], &tens.data_[0], dims_ * rows_ * cols_ * sizeof(T));
    }
    // move constructor
    Tensor(Tensor&& tens) noexcept
        : dims_(tens.dims_),
          rows_(tens.rows_),
          cols_(tens.cols_),
          data_(std::move(tens.data_)) {}
    // destructor
    ~Tensor() = default;
    // assignment operator
    Tensor& operator=(const Tensor& tens) {
        if(&tens != this) {
            memcpy(&data_[0], &tens.data_[0],
                   dims_ * rows_ * cols_ * sizeof(T));
        }
        return *this;
    }
    // move assignment operator
    Tensor& operator=(Tensor&& tens) noexcept {
        dims_ = tens.dims_;
        rows_ = tens.rows_;
        cols_ = tens.cols_;
        data_ = std::move(tens.data_);
        return *this;
    }
    T operator()(std::size_t dims, std::size_t row, std::size_t col) const {
        return data_[dims * rows_ * cols_ + row * cols_ + col];
    }
    T& operator()(std::size_t dims, std::size_t row, std::size_t col) {
        return data_[dims * rows_ * cols_ + row * cols_ + col];
    }
    bool operator==(const Tensor& tens) const;

   private:
    std::size_t dims_, rows_, cols_;
    std::vector<T> data_;
};  // class Tensor

using float_t = float;
using Matrixf = Matrix<float_t>;
using Matrixi = Matrix<int>;
using Tensorf = Tensor<float_t>;

}  // namespace coati
#endif
