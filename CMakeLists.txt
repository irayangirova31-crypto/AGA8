cmake_minimum_required(VERSION 3.16)
project(aga8_pvt_calculator LANGUAGES CXX)

# Стандарт языка C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Генерация compile_commands.json для VS Code
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Статическая библиотека
add_library(aga8_pvt_calculator STATIC
    src/aga8_pvt_calculator.cpp
    include/aga8_tables.h
    src/aga8_tables.cpp
)

# Исполняемый файл для тестирования
add_executable(test_aga8_pvt
    main.cpp          
)


# Пути к заголовочным файлам
target_include_directories(aga8_pvt_calculator
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)
# Для тестового приложения - только include директории
target_include_directories(test_aga8_pvt
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Линкуем библиотеку к тестовому приложению
target_link_libraries(test_aga8_pvt
    aga8_pvt_calculator
)

# Флаги компилятора
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(aga8_pvt_calculator PRIVATE -Wall -Wextra -Wpedantic)
     target_compile_options(test_aga8_pvt PRIVATE -Wall -Wextra -Wpedantic)
endif()