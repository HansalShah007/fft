cmake_minimum_required(VERSION 3.14)
project(hw5 LANGUAGES CXX)

find_package(UPCXX REQUIRED)




# cmake_minimum_required(VERSION 3.12)  # Update this line to avoid deprecation warning
project(iterative_fft)

# Fetch pybind11 from its GitHub repository
include(FetchContent)
FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG        v2.9.2  # Optionally specify a tag or commit
)
FetchContent_MakeAvailable(pybind11)

pybind11_add_module(iterative_fft /global/homes/b/baurkath/CS5220/project/v1/iterative_fft.cpp)


# Build the fft_upcpp executable
add_executable(fft_upcpp fft_upcpp.cpp)
target_link_libraries(fft_upcpp PRIVATE UPCXX::upcxx)
target_compile_definitions(fft_upcpp PRIVATE "KMER_LEN=19")


#add_executable(main_fftw main_fftw.cpp)


# Copy the job scripts
configure_file(job_fft_upcpp job_fft_upcpp COPYONLY)