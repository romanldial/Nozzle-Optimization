cmake_minimum_required(VERSION 3.20)

project(nozzle
    VERSION 0.1.0
    DESCRIPTION "Preliminary rocket nozzle design and optimization tool"
    LANGUAGES CXX)

# --- The reusable physics library (header-only) ---
add_library(nozzle_core INTERFACE)
target_include_directories(nozzle_core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(nozzle_core INTERFACE cxx_std_20)

# --- The demo executable ---
add_executable(nozzle_demo src/main.cpp)
target_link_libraries(nozzle_demo PRIVATE nozzle_core)