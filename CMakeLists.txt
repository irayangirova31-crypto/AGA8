cmake_minimum_required(VERSION 3.16)
project(aga8_pvt_calculator LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# 1. Статическая библиотека
add_library(aga8_pvt_calculator STATIC
    src/aga8_pvt_calculator.cpp
    src/aga8_tables.cpp
)

target_include_directories(aga8_pvt_calculator
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# 2. Демонстрационное приложение (из main.cpp)
add_executable(aga8_demo main.cpp)
target_link_libraries(aga8_demo PRIVATE aga8_pvt_calculator)

# 3. Набор тестов (из test_aga8.cpp)
add_executable(aga8_tests test_aga8.cpp)
target_link_libraries(aga8_tests PRIVATE aga8_pvt_calculator)

# Флаги компилятора
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(aga8_pvt_calculator PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(aga8_demo PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(aga8_tests PRIVATE -Wall -Wextra -Wpedantic)
endif()