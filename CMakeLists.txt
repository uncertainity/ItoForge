cmake_minimum_required(VERSION 3.16)
project(QuantEngine CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# This helps VS Code IntelliSense a lot
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(QuantLib REQUIRED)

# ---- Library (your engine code) ----
add_library(qe
    src/params.cpp
    src/path.cpp
    src/pricing.cpp
    src/surface.cpp
)

target_include_directories(qe
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(qe
    PUBLIC QuantLib::QuantLib
)

# ---- Executable (your main) ----
add_executable(demo
    app/main.cpp
)

target_link_libraries(demo
    PRIVATE qe
)