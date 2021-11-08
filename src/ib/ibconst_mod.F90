MODULE ibconst_mod
    USE core_mod, ONLY: realk, intk

    IMPLICIT NONE (type, external)

    ! TODO: consider moving these to ibcore_mod and read them from
    ! parameters.json like openaccur
    REAL(realk), PARAMETER :: maccur = 1.0e-4
    INTEGER(intk), PARAMETER :: nloopmax = 300
    INTEGER(intk), PARAMETER :: itermax = 10000
END MODULE ibconst_mod
