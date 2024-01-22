#pragma once
#define _USE_MATH_DEFINES
#include <complex>
#include <cmath>
#include <vector>

/// @brief Discrete Fourier Transform (by definition)
std::complex<double> DFT(const std::vector<std::complex<double>> &f, size_t k);

/// @brief Discrete Fourier Transform (by definition) but for all elements. Made for fair performance comparison with FFT.
std::vector<std::complex<double>> DFT(const std::vector<std::complex<double>> &f);

/// @brief Fast Fourier Transform (Cooley-Tukey radix-2)
/// @param f input data vector. Size of the vector must be a power of 2.
std::vector<std::complex<double>> FFT(const std::vector<std::complex<double>> &f);
