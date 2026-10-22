include_guard(GLOBAL)

set(FEDDLib_TEST_RUNTIME High CACHE STRING
    "Maximum expected test runtime tier: Low (<10 s), Medium (<30 s), High (all)")
set_property(CACHE FEDDLib_TEST_RUNTIME PROPERTY STRINGS Low Medium High)
set(_FEDDLib_RUNTIME_TIERS Low Medium High)
if(NOT FEDDLib_TEST_RUNTIME IN_LIST _FEDDLib_RUNTIME_TIERS)
    message(FATAL_ERROR
        "Invalid FEDDLib_TEST_RUNTIME='${FEDDLib_TEST_RUNTIME}'. Choose Low, Medium, or High.")
endif()

# Call after registering all tests in a directory. DEFAULT applies to that
# directory's tests; LOW/MEDIUM/HIGH override individual full CTest names.
# These are fixed laptop-based estimates, including every phase/MPI launch of
# a CTest case. They do not impose timeouts or read machine-local timing history.
function(feddlib_set_test_runtimes)
    cmake_parse_arguments(runtime "" "DEFAULT" "LOW;MEDIUM;HIGH" ${ARGN})
    if(runtime_UNPARSED_ARGUMENTS OR runtime_KEYWORDS_MISSING_VALUES)
        message(FATAL_ERROR "Invalid arguments to feddlib_set_test_runtimes in ${CMAKE_CURRENT_SOURCE_DIR}")
    endif()
    if(NOT runtime_DEFAULT)
        set(runtime_DEFAULT High)
    endif()
    set(tiers Low Medium High)
    if(NOT runtime_DEFAULT IN_LIST tiers)
        message(FATAL_ERROR "Invalid default test runtime tier '${runtime_DEFAULT}'")
    endif()
    list(FIND tiers "${FEDDLib_TEST_RUNTIME}" maximum_tier)
    get_property(tests DIRECTORY PROPERTY TESTS)
    foreach(test IN LISTS tests)
        set(tier "${runtime_DEFAULT}")
        set(overrides "")
        foreach(candidate IN LISTS tiers)
            string(TOUPPER "${candidate}" keyword)
            if(test IN_LIST runtime_${keyword})
                list(APPEND overrides "${candidate}")
                set(tier "${candidate}")
            endif()
        endforeach()
        list(LENGTH overrides override_count)
        if(override_count GREATER 1)
            message(FATAL_ERROR "Multiple runtime tiers assigned to test '${test}'")
        endif()

        # Append so package, RestartTests and Performance labels are retained.
        set_property(TEST "${test}" APPEND PROPERTY LABELS "Runtime${tier}")
        set_property(GLOBAL APPEND PROPERTY FEDDLib_RUNTIME_CLASSIFIED_TESTS "${test}")
        list(FIND tiers "${tier}" test_tier)
        if(test_tier GREATER maximum_tier)
            set_property(TEST "${test}" PROPERTY DISABLED TRUE)
        endif()
    endforeach()
endfunction()

# Require new test directories/cases to be classified too. Inspection can be
# recursive, but test properties must be set in their registration directory
# with the project's supported CMake versions.
function(feddlib_check_test_runtimes directory)
    get_property(tests DIRECTORY "${directory}" PROPERTY TESTS)
    get_property(classified GLOBAL PROPERTY FEDDLib_RUNTIME_CLASSIFIED_TESTS)
    foreach(test IN LISTS tests)
        if(NOT test IN_LIST classified)
            message(FATAL_ERROR
                "Test '${test}' has no runtime tier. Call feddlib_set_test_runtimes() after registering tests in ${directory}/CMakeLists.txt.")
        endif()
    endforeach()
    get_property(children DIRECTORY "${directory}" PROPERTY SUBDIRECTORIES)
    foreach(child IN LISTS children)
        feddlib_check_test_runtimes("${child}")
    endforeach()
endfunction()
