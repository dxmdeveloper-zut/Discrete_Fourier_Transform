#include "DFT.hpp"

#define SC_DBL(x) static_cast<double>(x)
inline
std::complex<double> DFT(const std::vector<std::complex<double>> &f, size_t k)
{
    const auto i = std::complex<double>(0, 1);
    const auto N = f.size();
    std::complex<double> ck = 0;

    for (int n = 0; n < N; n++) {
        ck += f[n] * std::exp(-2 * M_PI * i * SC_DBL(k) * SC_DBL(n) / SC_DBL(N));
    }
    return ck;
}

std::vector<std::complex<double>> FFT(const std::vector<std::complex<double>> &f)
{
    size_t N = f.size();
    if (N == 0) return {};
    if (N == 1) return f;

    std::vector<std::complex<double>> even(N / 2);
    std::vector<std::complex<double>> odd(N / 2);
    std::vector<std::complex<double>> FT(N);

    for (size_t i = 0; i < N / 2; i++) {
        even[i] = f[2 * i];
        odd[i] = f[2 * i + 1];
    }

    even = FFT(even);
    odd = FFT(odd);

    for (size_t i = 0; i < N / 2; i++) {
        auto zeta = std::polar(1.0, -2.0 * M_PI * SC_DBL(i) / SC_DBL(N)); // e^(-2πi/N)
        FT[i] = even[i] + zeta * odd[i];
        FT[i + N / 2] = even[i] - zeta * odd[i];
    }

    return FT;
}

std::vector<std::complex<double>> DFT(const std::vector<std::complex<double>> &f)
{
    size_t N = f.size();
    if (N == 0) return {};
    if (N == 1) return f;

    std::vector<std::complex<double>> FT(f.size());

    for (size_t k = 0; k < f.size(); k++) {
        FT[k] = DFT(f, k);
    }

    return FT;
}