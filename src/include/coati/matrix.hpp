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

class Matrix {
   public:
    Matrix(unsigned rows, unsigned cols, float value = 0.0f);
    Matrix(unsigned rows, unsigned cols, Matrix64f& eigen_m);
    Matrix(const Matrix& mat);                 // copy constructor
    Matrix(Matrix&& mat) noexcept;             // move constructor
    ~Matrix() = default;                       // destructor
    Matrix& operator=(const Matrix& mat);      // assignment operator
    Matrix& operator=(Matrix&& mat) noexcept;  // move assignment operator
    float operator()(unsigned row, unsigned col) const;
    float& operator()(unsigned row, unsigned col);
    bool operator==(const Matrix& mat) const;

   private:
    unsigned rows_, cols_;
    std::vector<float> data_;
};

inline Matrix::Matrix(unsigned rows, unsigned cols, float value)
    : rows_(rows), cols_(cols) {
    data_.resize(rows_ * cols_);
    for(unsigned i = 0; i < rows_; i++) {
        for(unsigned j = 0; j < cols_; j++) {
            data_[i * cols_ + j] = value;
        }
    }
}

inline Matrix::Matrix(unsigned rows, unsigned cols, Matrix64f& eigen_m)
    : rows_(rows), cols_(cols) {
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m(
        eigen_m);
    data_.resize(rows_ * cols_);
    for(unsigned i = 0; i < rows_ * cols_; i++) {
        data_[i] = m(i);
    }
}

inline Matrix::Matrix(const Matrix& mat) : rows_(mat.rows_), cols_(mat.cols_) {
    memcpy(&data_[0], &mat.data_[0], rows_ * cols_ * sizeof(float));
}

inline Matrix::Matrix(Matrix&& mat) noexcept
    : rows_(mat.rows_), cols_(mat.cols_), data_(std::move(mat.data_)){};

inline Matrix& Matrix::operator=(const Matrix& mat) {
    if(&mat != this) {
        memcpy(&data_[0], &mat.data_[0], rows_ * cols_ * sizeof(float));
    }
    return *this;
}

inline Matrix& Matrix::operator=(Matrix&& mat) noexcept {
    rows_ = mat.rows_;
    cols_ = mat.cols_;
    data_ = std::move(mat.data_);
    return *this;
};

inline float Matrix::operator()(unsigned row, unsigned col) const {
    return data_[row * cols_ + col];
}

inline float& Matrix::operator()(unsigned row, unsigned col) {
    return data_[row * cols_ + col];
}

inline bool Matrix::operator==(const Matrix& mat) const {
    if((rows_ != mat.rows_) || (cols_ != mat.cols_)) {
        return false;
    }

    for(unsigned i = 0; i < rows_ * cols_; i++) {
        if(data_[i] != mat.data_[i]) {
            return false;
        }
    }
    return true;
}

////////////////////////////////////////////////////////////////////////////////

class Tensor {
   public:
    Tensor(unsigned dims, unsigned rows, unsigned cols, float value = 0.0f);
    Tensor(const Tensor& tens);                 // copy constructor
    Tensor(Tensor&& tens) noexcept;             // move constructor
    ~Tensor() = default;                        // destructor
    Tensor& operator=(const Tensor& tens);      // assignment operator
    Tensor& operator=(Tensor&& tens) noexcept;  // move assignment operator
    float operator()(unsigned dims, unsigned rows, unsigned cols) const;
    float& operator()(unsigned dims, unsigned rows, unsigned cols);
    bool operator==(const Tensor& tens) const;

   private:
    unsigned dims_, rows_, cols_;
    std::vector<float> data_;
};

inline Tensor::Tensor(unsigned dims, unsigned rows, unsigned cols, float value)
    : dims_(dims), rows_(rows), cols_(cols) {
    data_.resize(dims_ * rows_ * cols_);
    for(unsigned i = 0; i < dims_ * rows_ * cols_; i++) {
        data_[i] = value;
    }
}

inline Tensor::Tensor(const Tensor& tens)
    : dims_(tens.dims_), rows_(tens.rows_), cols_(tens.cols_) {
    memcpy(&data_[0], &tens.data_[0], dims_ * rows_ * cols_ * sizeof(float));
}

inline Tensor::Tensor(Tensor&& tens) noexcept
    : dims_(tens.dims_),
      rows_(tens.rows_),
      cols_(tens.cols_),
      data_(std::move(tens.data_)) {}

inline Tensor& Tensor::operator=(const Tensor& tens) {
    if(&tens != this) {
        memcpy(&data_[0], &tens.data_[0],
               dims_ * rows_ * cols_ * sizeof(float));
    }
    return *this;
}

inline Tensor& Tensor::operator=(Tensor&& tens) noexcept {
    dims_ = tens.dims_;
    rows_ = tens.rows_;
    cols_ = tens.cols_;
    data_ = std::move(tens.data_);
    return *this;
}

inline float& Tensor::operator()(unsigned dim, unsigned row, unsigned col) {
    return data_[dim * rows_ * cols_ + row * cols_ + col];
}

inline float Tensor::operator()(unsigned dim, unsigned row,
                                unsigned col) const {
    return data_[dim * rows_ * cols_ + row * cols_ + col];
}
#endif
