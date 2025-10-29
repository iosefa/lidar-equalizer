cmake_minimum_required(VERSION 3.15)

project(lidar_equalizer VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall")

find_package(PDAL REQUIRED)

add_executable(lidar_equalizer
    src/main.cpp
    src/Equalizer.cpp
    src/OverlapFlagger.cpp
    src/OverlapOps.cpp
    src/OverlapSelfTest.cpp
)

target_include_directories(lidar_equalizer
    PRIVATE
        ${PDAL_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(lidar_equalizer
    PRIVATE
        ${PDAL_LIBRARIES}
)

include(CTest)
if (BUILD_TESTING)
    set(_test_input ${CMAKE_SOURCE_DIR}/data/5QKC570100.copc.laz)
    if (NOT EXISTS ${_test_input})
        set(_test_input ${CMAKE_SOURCE_DIR}/data/Lidar2018.copc.laz)
    endif()

    if (NOT EXISTS ${_test_input})
        message(WARNING "No COPC test fixture found; skipping CLI tests.")
    else()
        add_test(
            NAME lidar_equalizer_cli
            COMMAND
                ${CMAKE_COMMAND}
                    -DLIDAR_EQUALIZER_EXE=$<TARGET_FILE:lidar_equalizer>
                    -DINPUT_FILE=${_test_input}
                    -DOUTPUT_FILE=${CMAKE_BINARY_DIR}/tests/5QKC570100_equalized.laz
                    -P ${CMAKE_SOURCE_DIR}/tests/run_cli_test.cmake
        )

        add_test(
            NAME lidar_equalizer_cli_nonground
            COMMAND
                ${CMAKE_COMMAND}
                    -DLIDAR_EQUALIZER_EXE=$<TARGET_FILE:lidar_equalizer>
                    -DINPUT_FILE=${_test_input}
                    -DOUTPUT_FILE=${CMAKE_BINARY_DIR}/tests/5QKC570100_equalized_nonground.laz
                    -DCLASS_SCOPE=nonground
                    -P ${CMAKE_SOURCE_DIR}/tests/run_cli_test.cmake
        )
        add_test(
            NAME lidar_equalizer_cli_overlap
            COMMAND
                ${CMAKE_COMMAND}
                    -DLIDAR_EQUALIZER_EXE=$<TARGET_FILE:lidar_equalizer>
                    -DINPUT_FILE=${_test_input}
                    -DOUTPUT_FILE=${CMAKE_BINARY_DIR}/tests/5QKC570100_equalized_overlap.laz
                    -DFLAG_OVERLAP_ONLY=ON
                    -DMIN_POINTS_PER_PSID=5
                    -DOVERLAP_DILATE=1
                    -DSWATH_KEY=psid-channel
                    -P ${CMAKE_SOURCE_DIR}/tests/run_cli_test.cmake
        )
        add_test(
            NAME overlap_self_test
            COMMAND $<TARGET_FILE:lidar_equalizer> --self-test-overlap
        )
    endif()
endif()
