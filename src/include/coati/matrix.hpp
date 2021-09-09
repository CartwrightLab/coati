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

#include <Eigen/Dense>
#include <vector>

using Matrix64f = Eigen::Matrix<float, 64, 64>;

namespace coati {

using float_t = float;

template <class T>
class Matrix {
   public:
    Matrix(std::size_t rows, std::size_t cols, T value = 0.0f)
        : rows_(rows), cols_(cols) {
        data_.resize(rows_ * cols_);
        data_.assign(data_.size(), value);
    }
    Matrix(std::size_t rows, std::size_t cols, Matrix64f& eigen_m)
        : rows_(rows), cols_(cols) {
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m(
            eigen_m);
        data_.resize(rows_ * cols_);
        for(std::size_t i = 0; i < rows_ * cols_; i++) {
            data_[i] = m(i);
        }
    }
    // copy constructor
    Matrix(const Matrix& mat) : rows_(mat.rows_), cols_(mat.cols_) {
        memcpy(&data_[0], &mat.data_[0], rows_ * cols_ * sizeof(T));
    }
    // move constructor
    Matrix(Matrix&& mat) noexcept
        : rows_(mat.rows_), cols_(mat.cols_), data_(std::move(mat.data_)){};
    // destructor
    ~Matrix() = default;
    // assignment operator
    Matrix& operator=(const Matrix& mat) {
        if(&mat != this) {
            memcpy(&data_[0], &mat.data_[0], rows_ * cols_ * sizeof(T));
        }
        return *this;
    }
    // move assignment operator
    Matrix& operator=(Matrix&& mat) noexcept {
        rows_ = mat.rows_;
        cols_ = mat.cols_;
        data_ = std::move(mat.data_);
        return *this;
    }
    T operator()(std::size_t row, std::size_t col) const {
        return data_[row * cols_ + col];
    }
    T& operator()(std::size_t row, std::size_t col) {
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

   private:
    std::size_t rows_, cols_;
    std::vector<T> data_;
};  // class matrix

////////////////////////////////////////////////////////////////////////////////

// template <class T>
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

using Matrixf = Matrix<float_t>;
using Matrixi = Matrix<int>;
using Tensorf = Tensor<float_t>;

}  // namespace coati
#endif
