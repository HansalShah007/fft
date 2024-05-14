 

setup for the first time:

module load conda 
conda create --name parallel python=3.8
conda activate parallel
pip install -r requirements.txt 
 


mkdir build 
cd build


modify iterative_fft directory in CMakeLists.txt
move test.py into build + modify audio files directory 

cmake ..
make 


python ./test.py




module load contrib upcxx










run again:

module load conda 
module load contrib upcxx

cd build 
make 


sbatch job_fft_upcpp


python ./test.py






example for running file that includes fftw - modify path:


module load spack
spack list fftw
spack install fftw@3.3.10


g++ main_fftw.cpp -o main_fftw -I/global/homes/b/baurkath/spack-workspace/perlmutter/software/linux-sles15-zen3/gcc-12.3.0/fftw-3.3.10-xh4mphawnx4ag64b5w675rmxrz3pwxux/include -L/global/homes/b/baurkath/spack-workspace/perlmutter/software/linux-sles15-zen3/gcc-12.3.0/fftw-3.3.10-xh4mphawnx4ag64b5w675rmxrz3pwxux/lib -lfftw3


modify  job_fft_upcpp to have srun ./main_fftw
sbatch job_fft_upcpp









