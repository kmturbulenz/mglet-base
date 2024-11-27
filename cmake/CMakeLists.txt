include(CheckCSourceRuns)
include(CheckFortranSourceRuns)

include(Subdirlist.cmake)
include(TestQsort.cmake)

# testing Fortran language features
include(TestFortran2018Capability.cmake)
include(TestAbstractInterface.cmake)
include(TestElemental.cmake)

# testing C inter-operability
include(TestBindingC.cmake)
include(TestDataPtrC.cmake)
include(TestFunctionPtrC.cmake)

include(TestSelectRank.cmake)

# Test C++ headers
include(CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX("charconv" HAVE_CHARCONV)
if(NOT HAVE_CHARCONV)
    message(WARNING "C++17 include charconv not available")
endif()

