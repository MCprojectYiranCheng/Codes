cmake_minimum_required(VERSION 3.0)
project(MC)
set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

#find_package(GSL REQUIRED)



set(MC_SOURCES "${CMAKE_SOURCE_DIR}/src")
set(MC_INCLUDES "${CMAKE_SOURCE_DIR}/include")
set(MC_DOC "${CMAKE_SOURCE_DIR}/doc") 
set(MC_BIN "${CMAKE_SOURCE_DIR}/bin")
set(EXECUTABLE_OUTPUT_PATH ${MC_BIN})
#set(LIBS_GSL ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})


include_directories(${MC_INCLUDES} ${CMAKE_SYSTEM_INCLUDE_PATH})
#include_directories(${GSL_INCLUDE_DIRS} ${GSLCBLAS_INCLUDE_DIRS})

add_executable(asianOptionNew "${MC_SOURCES}/asianOptionNew.cpp" "${MC_SOURCES}/asa241.cpp" "${MC_SOURCES}/importanceSamplingAndVectorTools.cpp" "${MC_SOURCES}/stratifiedMCClass.cpp")
add_executable(simpleNormal "${MC_SOURCES}/simpleNormal.cpp" "${MC_SOURCES}/asa241.cpp" "${MC_SOURCES}/importanceSamplingAndVectorTools.cpp" "${MC_SOURCES}/stratifiedMCClass.cpp")


#target_link_libraries(test_beta ${LIBS_GSL})

