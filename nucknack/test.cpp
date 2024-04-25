#include <iostream>
#include <stdio.h>
#include <omp.h>
#include <chrono>
#include <thread>
#include <vector>
#include <complex>

using namespace std;

int reverse(int num, int lg_n) {
	int res = 0;
	for (int i = 0; i < lg_n; i++) {
			if (num & (1 << i))
					res |= 1 << (lg_n - 1 - i);
	}
	return res;
}

int main() {
	cout << reverse(5, 4) << endl;

	int lg_n = 3;
	for (int i = 0; i < 8; i++) {
		int reversed_bits = reverse(i, lg_n);
		cout << i << ' ' << reversed_bits << endl;
			if (i < reversed_bits)
				cout << ' ' << i << ' ' << reversed_bits << endl;
				//	swap(a[i], a[reversed_bits]);
		cout << endl;
	}

	return 0;
}