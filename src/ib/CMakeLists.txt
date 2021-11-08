add_library(ib STATIC
    blockbp_mod.F90
    blockcheck_mod.F90
    blockclosetoboundary_mod.F90
    blockluecken_mod.F90
    blocknodes_mod.F90
    bubvbw_mod.F90
    calcauavaw_mod.F90
    calcfacearea_mod.F90
    calcnormals_mod.F90
    checkzelle_mod.F90
    ctof_mod.F90
    cutcorner_mod.F90
    fillfollower_mod.F90
    filling_mod.F90
    findinterface_mod.F90
    flzelle_mod.F90
    freekante_mod.F90
    freepressure_mod.F90
    ftoc_mod.F90
    gc_blockbpfeld_mod.F90
    gc_blockbp_mod.F90
    gc_blockface_mod.F90
    gc_createstencils_mod.F90
    gc_finishknotenbzelltyp_mod.F90
    gc_mod.F90
    gc_restrict_mod.F90
    gc_stencils_mod.F90
    gc_totwasser_mod.F90
    gc_zelltyp_mod.F90
    ib_mod.F90
    ibconst_mod.F90
    ibcore_mod.F90
    ibmodel_mod.F90
    knotenundkanten_mod.F90
    noib_mod.F90
    noib_restrict_mod.F90
    openbubvbw_mod.F90
    parent_mod.F90
    punktekoordinaten_mod.F90
    stencils_mod.F90
    topol_mod.F90
)

target_include_directories(ib
    PRIVATE ${MPI_Fortran_INCLUDE_PATH}
    PRIVATE ${HDF5_Fortran_INCLUDE_DIRS}
)

target_link_libraries(ib
    PRIVATE core
    PRIVATE ${MPI_Fortran_LIBRARIES}
)

set_source_files_properties(checkzelle_mod.F90 PROPERTIES COMPILE_FLAGS
    $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-Wno-compare-reals>
)

set_source_files_properties(cutcorner_mod.F90 PROPERTIES COMPILE_FLAGS
    $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-Wno-compare-reals>
)

set_source_files_properties(ftoc_mod.F90 PROPERTIES COMPILE_FLAGS
    $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-Wno-compare-reals>
)
