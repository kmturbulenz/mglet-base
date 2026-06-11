MODULE gpuruntime_mod
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_char, c_null_char, c_int

    IMPLICIT NONE (type, external)
    PRIVATE

    INTERFACE
        SUBROUTINE gpu_runtime_init_c() BIND(C, name="gpu_backend_runtime_init")
            ! Subroutine arguments
            ! none...
        END SUBROUTINE gpu_runtime_init_c
        
        SUBROUTINE gpu_runtime_finalize_c() &
            BIND(C, name="gpu_backend_runtime_finalize")
            ! Subroutine arguments
            ! none...
        END SUBROUTINE gpu_runtime_finalize_c
    END INTERFACE

    PUBLIC :: gpu_runtime_init, gpu_runtime_finalize
CONTAINS
    SUBROUTINE gpu_runtime_init()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        CALL gpu_runtime_init_c()
    END SUBROUTINE gpu_runtime_init


    SUBROUTINE gpu_runtime_finalize()
        ! Subroutine arguments
        ! none...

        ! Local variables
        ! none...

        CALL gpu_runtime_finalize_c()
    END SUBROUTINE gpu_runtime_finalize
END MODULE gpuruntime_mod
