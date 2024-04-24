#include <iostream>
#include <stdio.h>
#include <omp.h>
#include <chrono>
#include <thread>

using namespace std;

//  g++ -o omp_fft fft_omp.cpp -fopenmp
/*
export OMP_NUM_THREADS=32
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
*/

int main() {
	auto start = std::chrono::steady_clock::now();
	#pragma omp master
	{
		std::cout << "Threads: " << omp_get_num_threads() << endl;
	}

    int shared_variable = 0;
    #pragma omp parallel shared(shared_variable)
    {
        // Increment the shared variable with a lock
        #pragma omp atomic
        shared_variable++;

		std::this_thread::sleep_for(std::chrono::seconds(1));
    }

	#pragma omp barrier

	#pragma omp master
	{
		cout << "Shared " << shared_variable << endl;
		auto end = std::chrono::steady_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
		std::cout << "Elapsed time: " << 1.0*duration/1000000 << " sec" << std::endl;
	}
    return 0;
}
