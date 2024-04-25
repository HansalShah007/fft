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
    int shared_variable = 0;
	#pragma omp master
	{
		#pragma omp parallel shared(shared_variable)
		{
			cout << omp_get_thread_num() << endl;
			//shared_variable += omp_get_thread_num();
		}
	}
	
	/*
	#pragma omp barrier

	#pragma omp master
	{
		cout << shared_variable << endl;
	}
	*/
    return 0;
}
