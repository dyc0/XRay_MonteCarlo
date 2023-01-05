cmake_minimum_required(VERSION 3.10)

project(XRAY_SIM VERSION 1.0)


set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

option(DEBUGGING "Toggle debugging mode for symbol generation etc." ON)
option(ENABLE_TESTING "Choose whether to compile test files." ON)

if(DEBUGGING)
  set(CMAKE_BUILD_TYPE Debug)
endif()

add_library(compiler_flags INTERFACE)
target_compile_features(compiler_flags INTERFACE cxx_std_11)

if(DEBUGGING)
  target_compile_options(compiler_flags INTERFACE
    -Wall -Wextra -Wshadow -Wformat=2 -Wunused -mavx
  )
else()
  target_compile_options(compiler_flags INTERFACE
    -O0 -mavx2
  )
endif()

get_target_property(CompFl compiler_flags CMAKE_C_FLAGS)
message(${CompFl})

add_subdirectory(LIBRARY)
list(APPEND XR_LIBS XRaySimulator)

if(ENABLE_TESTING)
  add_subdirectory(TESTS)
endif()

#add_executable(solution solution.cpp)
#target_link_libraries(solution ${XR_LIBS} compiler_flags)
