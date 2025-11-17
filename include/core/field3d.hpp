#pragma once

#include <memory>
#include <stdexcept>

namespace fvm3d::core {

/**
 * 3D field storage with Structure-of-Arrays (SoA) memory layout.
 * Data is stored as: [var][i][j][k]
 * This layout improves cache efficiency for SIMD operations.
 */
template<typename T = double>
class Field3D {
public:
    /**
     * Constructor: allocate field with specified dimensions.
     * @param nvars: number of variables per cell (e.g., 5 for Euler: rho, rho_u, rho_v, rho_w, E)
     * @param nx, ny, nz: grid dimensions (total including ghost cells)
     */
    Field3D(int nvars, int nx, int ny, int nz)
        : nvars_(nvars), nx_(nx), ny_(ny), nz_(nz) {
        if (nvars <= 0 || nx <= 0 || ny <= 0 || nz <= 0) {
            throw std::invalid_argument("Field3D dimensions must be positive");
        }
        size_ = (size_t)nvars * nx * ny * nz;
        data_ = std::make_unique<T[]>(size_);
    }

    // Move semantics (copy is deleted)
    Field3D(Field3D&& other) noexcept = default;
    Field3D& operator=(Field3D&& other) noexcept = default;

    // Delete copy
    Field3D(const Field3D&) = delete;
    Field3D& operator=(const Field3D&) = delete;

    // Accessors
    int nvars() const { return nvars_; }
    int nx() const { return nx_; }
    int ny() const { return ny_; }
    int nz() const { return nz_; }
    size_t size() const { return size_; }
    T* data() { return data_.get(); }
    const T* data() const { return data_.get(); }

    /**
     * Element access: field(variable, i, j, k)
     */
    inline T& operator()(int v, int i, int j, int k) {
        return data_[index(v, i, j, k)];
    }

    inline const T& operator()(int v, int i, int j, int k) const {
        return data_[index(v, i, j, k)];
    }

    /**
     * Fill entire field with a value.
     */
    void fill(T value) {
        for (size_t i = 0; i < size_; i++) {
            data_[i] = value;
        }
    }

    /**
     * Copy data from another field (requires same dimensions).
     */
    void copy_from(const Field3D& other) {
        if (nvars_ != other.nvars_ || nx_ != other.nx_ ||
            ny_ != other.ny_ || nz_ != other.nz_) {
            throw std::invalid_argument("Field dimensions mismatch");
        }
        std::copy(other.data_.get(), other.data_.get() + size_, data_.get());
    }

    /**
     * Assign data from another field element-wise (for temporal integration).
     * This allows *this = other; by copying data.
     */
    Field3D& assign(const Field3D& other) {
        copy_from(other);
        return *this;
    }

    // ========== Vectorized Operations for Time Integration ==========

    /**
     * AXPY operation: this = a * x + y
     * Vectorized operation on entire field (single pass through memory).
     */
    void axpy(T a, const Field3D& x, const Field3D& y) {
        #pragma omp simd
        for (size_t i = 0; i < size_; i++) {
            data_[i] = a * x.data_[i] + y.data_[i];
        }
    }

    /**
     * Linear combination: this = a * x + b * y
     * Vectorized operation on entire field.
     */
    void linear_combination(T a, const Field3D& x, T b, const Field3D& y) {
        #pragma omp simd
        for (size_t i = 0; i < size_; i++) {
            data_[i] = a * x.data_[i] + b * y.data_[i];
        }
    }

    /**
     * Three-term linear combination: this = a * x + b * y + c * z
     * Vectorized operation on entire field.
     */
    void linear_combination_3(
        T a, const Field3D& x,
        T b, const Field3D& y,
        T c, const Field3D& z
    ) {
        #pragma omp simd
        for (size_t i = 0; i < size_; i++) {
            data_[i] = a * x.data_[i] + b * y.data_[i] + c * z.data_[i];
        }
    }

    /**
     * Scaled add: this = this + a * x
     * In-place AXPY operation (BLAS DAXPY pattern).
     */
    void add_scaled(T a, const Field3D& x) {
        #pragma omp simd
        for (size_t i = 0; i < size_; i++) {
            data_[i] += a * x.data_[i];
        }
    }

    /**
     * Scale entire field: this = a * this
     */
    void scale(T a) {
        #pragma omp simd
        for (size_t i = 0; i < size_; i++) {
            data_[i] *= a;
        }
    }

private:
    int nvars_;
    int nx_, ny_, nz_;
    size_t size_;
    std::unique_ptr<T[]> data_;

    /**
     * Compute linear index for SoA layout.
     * Memory layout: [var][i][j][k] in row-major order
     */
    inline size_t index(int v, int i, int j, int k) const {
        return (size_t)v * (nx_ * ny_ * nz_) +
               (size_t)i * (ny_ * nz_) +
               (size_t)j * nz_ +
               (size_t)k;
    }
};

// Convenience typedef
using StateField3D = Field3D<double>;

} // namespace fvm3d::core
