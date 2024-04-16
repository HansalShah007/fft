cmake_minimum_required(VERSION 3.12)  # Update this line to avoid deprecation warning
project(fft)

# Fetch pybind11 from its GitHub repository
include(FetchContent)
FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG        v2.9.2  # Optionally specify a tag or commit
)
FetchContent_MakeAvailable(pybind11)

pybind11_add_module(fft /workspace/fft/fft.cpp)
