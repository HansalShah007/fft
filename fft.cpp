#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <vector>
#include <complex>
#include <cmath>

using namespace std;

const double PI = acos(-1);

// Function to find the nearest power of 2
int nearest_pow_2(int n) {
    int power = 1;
    while (power < n) power <<= 1;
    return power;
}

// Recursive FFT function
vector<complex<double>> fft(vector<complex<double>>& a) {
    int n = a.size();

    if (n == 1)
        return vector<complex<double>>(1, a[0]);

    // Check if the size is a power of 2, if not pad with zeros
    int n2 = nearest_pow_2(n);
    if (n2 > n) {
        a.resize(n2, 0);
        n = n2;
    }

    vector<complex<double>> even(n / 2), odd(n / 2);
    for (int i = 0; i < n / 2; ++i) {
        even[i] = a[i * 2];
        odd[i] = a[i * 2 + 1];
    }

    auto y_even = fft(even);
    auto y_odd = fft(odd);

    vector<complex<double>> y(n);
    complex<double> omega_n = exp(complex<double>(0, 2 * PI / n));
    complex<double> omega = 1;

    for (int j = 0; j < n / 2; ++j) {
        y[j] = y_even[j] + omega * y_odd[j];
        y[j + n / 2] = y_even[j] - omega * y_odd[j];
        omega *= omega_n;
    }

    return y;
}

// Recursive IFFT function
vector<complex<double>> ifft(vector<complex<double>>& a, int depth) {
    int n = a.size();

    if (n == 1)
        return vector<complex<double>>(1, a[0]);

    int n2 = nearest_pow_2(n);
    if (n2 > n) {
        a.resize(n2, 0);
        n = n2;
    }

    vector<complex<double>> even(n / 2), odd(n / 2);
    for (int i = 0; i < n / 2; ++i) {
        even[i] = a[i * 2];
        odd[i] = a[i * 2 + 1];
    }

    auto y_even = ifft(even, depth+1);
    auto y_odd = ifft(odd,depth+1);

    vector<complex<double>> y(n);
    complex<double> omega_n = exp(complex<double>(0, -2 * PI / n));
    complex<double> omega = 1;

    for (int j = 0; j < n / 2; ++j) {
        y[j] = y_even[j] + omega * y_odd[j];
        y[j + n / 2] = y_even[j] - omega * y_odd[j];
        omega *= omega_n;
    }

    if (depth==1){
        for (auto& value : y) value /= y.size();
    }

    return y;
}

PYBIND11_MODULE(fft, m) {
    m.doc() = "FFT module implemented in C++";
    m.def("fft", &fft, "A function which computes the FFT of a complex number vector");
    m.def("ifft", &ifft, "A function which computes the inverse FFT of a complex number vector");
    m.def("nearest_pow_2", &nearest_pow_2, "A function that returns the next nearest value that is a power of 2");
}