cmake_minimum_required(VERSION 3.0)
project(RQMC)
set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

find_package(GSL REQUIRED)



set(RQMC_SOURCES "${CMAKE_SOURCE_DIR}/src")
set(RQMC_INCLUDES "${CMAKE_SOURCE_DIR}/include")
set(RQMC_DOC "${CMAKE_SOURCE_DIR}/doc") 
set(RQMC_BIN "${CMAKE_SOURCE_DIR}/bin")
set(EXECUTABLE_OUTPUT_PATH ${RQMC_BIN})
set(LIBS_GSL ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})


include_directories(${RQMC_INCLUDES} ${CMAKE_SYSTEM_INCLUDE_PATH})
include_directories(${GSL_INCLUDE_DIRS} ${GSLCBLAS_INCLUDE_DIRS})

add_executable(test_beta "${RQMC_SOURCES}/main.cpp")


target_link_libraries(test_beta ${LIBS_GSL})

