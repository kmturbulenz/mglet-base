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
