cmake_minimum_required(VERSION 3.15)
set(PROJECT_NAME NSST)
project(${PROJECT_NAME})

set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_STANDARD 17)

set(HEADERS
        Config.h )

set(SOURCES
        src/NSST.cpp
        main.cpp )


include_directories(external_lib include)
add_executable(${PROJECT_NAME} ${HEADERS} ${SOURCES})

set(EXTERNAL_LIB_PATH
	${PROJECT_SOURCE_DIR}/external_lib)
	
set(EIGEN_PATH
	${EXTERNAL_LIB_PATH}/eigen-3.4.0)

message("-- Linking FFTW3")

if(WIN32)
    target_include_directories(${PROJECT_NAME} PRIVATE "${EXTERNAL_LIB_PATH}/fftw-3.3.5-dll64" "${EIGEN_PATH}")
    target_link_directories(${PROJECT_NAME} PRIVATE "${EXTERNAL_LIB_PATH}/fftw-3.3.5-dll64" "${EIGEN_PATH}")
    target_link_libraries(${PROJECT_NAME} PRIVATE libfftw3f-3)
endif(WIN32)

if(UNIX)
    target_include_directories(${PROJECT_NAME} PRIVATE "${EXTERNAL_LIB_PATH}/fftw-3.3.5" "${EIGEN_PATH}")
    target_link_directories(${PROJECT_NAME} PRIVATE "${EXTERNAL_LIB_PATH}/fftw-3.3.5" "${EIGEN_PATH}")
    target_link_libraries(${PROJECT_NAME} PRIVATE fftw3f)
endif(UNIX)
