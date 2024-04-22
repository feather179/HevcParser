cmake_minimum_required(VERSION 3.5)

project(HevcParser VERSION 0.1 LANGUAGES CXX)

set(EXE_NAME HevcParser)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB COMMON_SRC_FILES "common/*.cpp")
file(GLOB DECODER_SRC_FILES "decoder/*.cpp")
file(GLOB READER_SRC_FILES "reader/*.cpp")

set(SRC_FILES
    main.cpp
    ${COMMON_SRC_FILES}
    ${DECODER_SRC_FILES}
    ${READER_SRC_FILES}
)

add_executable(${EXE_NAME} ${SRC_FILES})

include_directories(common decoder reader)
