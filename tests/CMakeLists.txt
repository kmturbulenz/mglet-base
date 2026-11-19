enable_testing()

if(MGLET_OFFLOAD AND MGLET_OFFLOAD_BACKEND MATCHES "^(hip|cuda)$")
    set(sip_test_source sip_hyperplane_test.cu)
    if(MGLET_OFFLOAD_BACKEND STREQUAL "hip")
        set_source_files_properties(${sip_test_source} PROPERTIES LANGUAGE HIP)
    else()
        set_source_files_properties(${sip_test_source} PROPERTIES LANGUAGE CUDA)
    endif()
    add_executable(sip_hyperplane_test ${sip_test_source})
    set_target_properties(sip_hyperplane_test PROPERTIES
        CXX_STANDARD 17 HIP_STANDARD 17 CUDA_STANDARD 17
        LINKER_LANGUAGE Fortran)
    target_link_libraries(sip_hyperplane_test PRIVATE mglet_openmp
        ${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES}
        $<$<STREQUAL:${MGLET_OFFLOAD_BACKEND},hip>:hip::device>
        $<$<STREQUAL:${MGLET_OFFLOAD_BACKEND},cuda>:CUDA::cudart>)
    add_test(NAME sip_hyperplane COMMAND sip_hyperplane_test)
endif()

# By default more expensive testcases are also enabled - but can be disabled
# to make testing quicker
option(MGLET_TEST_Additional "Enable additional testcases" ON)
if(MGLET_TEST_Additional)
    message(STATUS "Enabling additional testcases")
    add_subdirectory("Additional-tests")
else()
    message(STATUS "Disabling additional testcases")
endif()

SUBDIRLIST(MGLET_TESTCASES ${CMAKE_CURRENT_SOURCE_DIR})
FOREACH(test ${MGLET_TESTCASES})
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${test}/run.sh)
        file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
        add_test(
            NAME ${test}
            COMMAND ${CMAKE_SOURCE_DIR}/ci/run-test.sh ${CMAKE_CURRENT_SOURCE_DIR}/${test} $<TARGET_FILE:mglet>
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test}
        )

        # Only working for CMake 3.27 and onwards - but earlier versions
        # does not seem to care
        set_tests_properties(${test} PROPERTIES TIMEOUT_SIGNAL_NAME SIGINT)
        set_tests_properties(${test} PROPERTIES TIMEOUT_SIGNAL_GRACE_PERIOD 4)
    endif()
ENDFOREACH()
