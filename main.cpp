#include <iostream>
#include <chrono>
#include "DFT.hpp"

//#define VERBOSE 1
#define PLOTTING 1
#ifdef PLOTTING
#include <matplotlibcpp.h>
#endif

#define HR_CLOCK_NOW() std::chrono::high_resolution_clock::now()
#define HR_CLOCK_DURATION std::chrono::duration_cast<std::chrono::milliseconds>

const size_t N = 1 << 14; // 2^14 = 16384

namespace plt = matplotlibcpp;

int main()
{
    // DFT test
    std::vector<std::complex<double>> samples(N);
    std::vector<std::complex<double>> fft_res(N);

    // Sample Generation
    const double freq = 23;
    const double sample_rate = 43.3;

    for (int i = 0; i < N; i++) {
        double t = i / sample_rate;
        //samples[i] = (4*i+24);
        samples[i] = std::sin(2.0 * M_PI * freq * t) * 1600;
#ifdef VERBOSE
        std::cout << "Sample[" << i << "] = " << samples[i] << std::endl;
#endif
    }

    // FFT, DFT benchmark
    auto fft_test_start = HR_CLOCK_NOW();
    fft_res = FFT(samples);
    auto fft_test_end = HR_CLOCK_NOW();

    auto dft_test_start = HR_CLOCK_NOW();
    auto dft_res = DFT(samples);
    auto dft_test_end = HR_CLOCK_NOW();

#ifdef VERBOSE
    for(size_t i = 0; i < N; i++) {
        std::cout << "FFT[" << i << "] = " << fft_res[i] << std::endl;
        std::cout << "DFT[" << i << "] = " << dft_res[i] << std::endl;
    }
#endif
    std::cout << "DFT for " << N << " samples took " << HR_CLOCK_DURATION(dft_test_end - dft_test_start).count()
              << " ms" << std::endl;
    std::cout << "FFT for " << N << " samples took " << HR_CLOCK_DURATION(fft_test_end - fft_test_start).count()
              << " ms" << std::endl;

    double err = 0.0;
    for (size_t i = 0; i < N; i++) {
        err += std::abs(fft_res[i] - dft_res[i]);
    }
    std::cout << "avg error: " << err / N << std::endl;

#ifdef PLOTTING
    // clear DFT result
    dft_res.clear();

    std::vector<double> dsamples(N);
    std::vector<double> ft_real(N), ft_img(N);
    std::vector<double> magnitude(N), phase(N);

    // Magnitude, Phase calculation
    for (size_t k = 0; k < N; k++) {
        dsamples[k] = samples[k].real();
        ft_real[k] = fft_res[k].real();
        ft_img[k] = fft_res[k].imag();
        magnitude[k] = std::sqrt(std::pow(ft_real[k], 2) + std::pow(ft_img[k], 2));
        phase[k] = std::atan2(ft_img[k], ft_real[k]);
    }

    // clear FFT results and samples
    fft_res.clear();
    samples.clear();

    // Plotting
    plt::plot(dsamples);
    plt::title("Samples");

    dsamples.clear();

    plt::figure();
    try { plt::subplot(1, 3, 1); } catch (...) {}
    plt::title("Amplitude (Fourier transform)");
    plt::plot(ft_real, {{"label", "real (cos)"}});
    plt::plot(ft_img, {{"label", "imaginary (sin)"}});
    plt::legend();

    ft_real.clear();
    ft_img.clear();

    try {
        plt::subplot(1, 3, 2);
    } catch (...) {
        plt::figure();
    }
    plt::title("Magnitude");
    plt::plot(magnitude);

    magnitude.clear();

    try {
        plt::subplot(1, 3, 3);
    } catch (...) {
        plt::figure();
    }
    plt::title("Phase");
    plt::plot(phase);

    phase.clear();

    plt::show();
#endif

    return 0;
}