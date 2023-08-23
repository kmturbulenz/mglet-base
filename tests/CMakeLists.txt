enable_testing()

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
