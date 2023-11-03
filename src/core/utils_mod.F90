MODULE utils_mod
    USE err_mod, ONLY: errr
    USE precision_mod, ONLY: intk, realk
    USE qsort_mod, ONLY: sortix
    USE comms_mod, ONLY: myid

    IMPLICIT NONE(type, external)
    PRIVATE

    INTERFACE
        ! void create_directory_f(const char *path)
        SUBROUTINE create_directory_f(path) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR
            CHARACTER(C_CHAR), INTENT(IN) :: path(*)
        END SUBROUTINE create_directory_f
    END INTERFACE

    PUBLIC :: unique, ind2sub, sub2ind, create_directory, &
        most_frequent_nonzero, read_datfile, get_stag_shift, get_idx

CONTAINS
    SUBROUTINE unique(result, input)
        ! Return an array with the unique elements of the input array

        ! Subroutine arguments
        INTEGER(intk), ALLOCATABLE, INTENT(out) :: result(:)
        INTEGER(intk), INTENT(in) :: input(:)

        ! Local variables
        INTEGER(intk) :: nuniq, ndata, i, this_elm, prev_elm
        INTEGER(intk), ALLOCATABLE :: idx(:)

        ! Useful in case user swaps arguments
        IF (ALLOCATED(result)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Number if input data elements
        ndata = SIZE(input)

        ! Shortcuts
        IF (ndata == 0) THEN
            ALLOCATE(result(0))
            RETURN
        ELSE IF (ndata == 1) THEN
            ALLOCATE(result(1))
            result(1) = input(1)
            RETURN
        END IF

        ALLOCATE(idx(ndata))

        CALL sortix(ndata, input, idx)

        ! There is at least always 1 unique element - the first one
        nuniq = 1
        DO i = 2, ndata
            this_elm = input(idx(i))
            prev_elm = input(idx(i-1))

            ! If the element was equal to the previous, it is not unique
            IF (this_elm /= prev_elm) THEN
                nuniq = nuniq + 1
            END IF
        END DO

        ALLOCATE(result(nuniq))
        result = 0
        result(1) = input(idx(1))

        nuniq = 1
        DO i = 2, ndata
            this_elm = input(idx(i))
            prev_elm = input(idx(i-1))

            ! If the element was equal to the previous, it is not unique
            IF (this_elm /= prev_elm) THEN
                nuniq = nuniq + 1
                result(nuniq) = this_elm
            END IF
        END DO
    END SUBROUTINE unique


    INTEGER(intk) FUNCTION most_frequent_nonzero(list)
        ! Return the most frequent notzero element in a list. If two elements
        ! have the same frequence of occurrence, return the biggest

        ! Function arguments
        INTEGER(intk), INTENT(in) :: list(:)

        ! Local variables
        INTEGER(intk) :: ndata, i, max_freq, curr_freq, curr_elm
        INTEGER(intk), ALLOCATABLE :: idx(:)

        ! Number if input data elements
        ndata = SIZE(list)

        ! Shortcuts
        IF (ndata == 0) THEN
            most_frequent_nonzero = 0
            RETURN
        ELSE IF (ndata == 1) THEN
            most_frequent_nonzero = list(1)
            RETURN
        END IF

        ALLOCATE(idx(ndata))
        CALL sortix(ndata, list, idx)

        most_frequent_nonzero = 0
        max_freq = 0
        i = 1
        DO WHILE (i <= ndata)
            curr_elm = list(idx(i))
            IF (curr_elm /= 0) THEN
                curr_freq = 1
            ELSE
                curr_freq = 0
            END IF

            DO WHILE (i+1 <= ndata)
                IF (list(idx(i+1)) == curr_elm) THEN
                    i = i + 1
                    IF (curr_elm /= 0) THEN
                        curr_freq = curr_freq + 1
                    END IF
                ELSE
                    EXIT
                ENDIF
            END DO

            ! Here always prefer the biggest if they are equal
            ! (hence <= instead of just < )
            IF (max_freq <= curr_freq) THEN
                max_freq = curr_freq
                most_frequent_nonzero = curr_elm
            END IF

            i = i + 1
        END DO

        DEALLOCATE(idx)
    END FUNCTION most_frequent_nonzero


    PURE ELEMENTAL SUBROUTINE ind2sub(ind, k, j, i, kk, jj, ii)
        ! Function converts a one dimensional index 'ind' into
        ! a tuple of three-dimensional indices (k, j, i) for an array.
        ! Size of the three-dimensional array must be known a priori.

        ! Subroutine arguments
        INTEGER(intk), INTENT(IN) :: ind
        INTEGER(intk), INTENT(OUT) :: i, j, k
        INTEGER(intk), INTENT(IN) :: ii, jj, kk

        k = MOD(ind-1, kk)+1
        j = MOD((ind-k)/kk, jj)+1
        i = (ind-k-(j-1)*kk)/(kk*jj)+1
    END SUBROUTINE ind2sub


    PURE ELEMENTAL SUBROUTINE sub2ind(ind, k, j, i, kk, jj, ii)
        ! Function converts a tuple of three-dimensional indices (k, j, i) into
        ! a one dimensional index \( ind \) for an array.
        ! Size of the three-dimensional array must be known a priori.

        ! Subroutine arguments
        INTEGER(intk), INTENT(out) :: ind
        INTEGER(intk), INTENT(in) :: i, j, k
        INTEGER(intk), INTENT(in) :: ii, jj, kk

        ind = k + (j-1)*kk + (i-1)*kk*jj
    END SUBROUTINE sub2ind


    SUBROUTINE create_directory(path)
        USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR, C_NULL_CHAR
        ! Check if a directory is present and create it if not

        ! Subroutine arguments
        CHARACTER(len=*), INTENT(in) :: path

        ! Local variables
        CHARACTER(C_CHAR), DIMENSION(LEN(path)+1) :: c_path

        ! Add trailing C_NULL_CHAR to key
        c_path = TRANSFER(path, c_path)
        c_path(LEN_TRIM(path)+1) = C_NULL_CHAR

        CALL create_directory_f(c_path)
    END SUBROUTINE create_directory


    SUBROUTINE read_datfile(array, filename)
        ! Reads a dat-file, i.e. a simple text file with data in columns. All
        ! data is assumed to be floating point and whitespace separated.
        !
        ! TODO: enhance functionality to allow '#' as a comment marker
        ! TODO: outsource to external library or C++ code?? The Fortran
        !       implementation here does not win any beauty award...

        ! Subroutine arguments
        REAL(realk), INTENT(out), ALLOCATABLE :: array(:, :)
        CHARACTER(len=*), INTENT(in) :: filename

        ! Local variables
        INTEGER(intk) :: nrow, ncol, i, j
        INTEGER(intk), PARAMETER :: nchar_max = 1024    ! Characters per line
        INTEGER(intk), PARAMETER :: ncol_max = 16       ! Max. number of columns
        INTEGER(intk), PARAMETER :: nrow_max = 1048576  ! Max. number of rows
        INTEGER :: unit, ierr
        LOGICAL :: fileexist

        CHARACTER(len=nchar_max) :: line
        REAL(realk) :: value

        ! Print sensible message if file does not exist
        INQUIRE(file=filename, EXIST=fileexist)
        IF (.NOT. fileexist) THEN
            WRITE(*,*) "File does not exist: ", filename
            CALL errr(__FILE__, __LINE__)
        END IF

        ! First step is to determine the number of columns. This is done by
        ! reading one line of data into a character string, then try to read
        ! values from that string until it fails, then you have determined the
        ! number of columns
        OPEN(NEWUNIT=unit, FILE=filename, ACTION='READ', IOSTAT=ierr)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        READ(unit, '(A)', iostat=ierr) line
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        ncol = 0
        DO i = 1, ncol_max + 1
            READ(line, *, iostat=ierr) (value, j=1, i)
            IF (ierr /= 0) THEN
                ncol = i - 1
                EXIT
            END IF
        END DO

        ! In case of more than ncol_max columns or just whitespace
        IF (ncol == 0) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! We have already read the first row from the file, let's continue
        ! reading to determine the number of rows
        DO i = 2, nrow_max + 1
            READ(unit, '(A)', iostat=ierr) line
            IF (ierr /= 0) THEN
                nrow = i - 1
                EXIT
            END IF
        END DO

        ! In case of more than nrow_max rows
        IF (nrow == 0) THEN
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Rewind to read from start, allocate array and read into it
        REWIND(unit)
        ALLOCATE(array(ncol, nrow))
        READ(unit, *, iostat=ierr) ((array(i, j), i=1,ncol), j=1,nrow)
        IF (ierr /= 0) CALL errr(__FILE__, __LINE__)

        CLOSE(unit)

    END SUBROUTINE read_datfile


    ! Compute shifts used for interpolation of different fields with different
    ! staggering properties
    ! istart, istop, must be initialized with grid shape
    ! (istart = 1, istop = ii)
    !
    ! stagout is the target staggering of the result
    ! stag1 is the input field staggering
    PURE SUBROUTINE get_stag_shift(i1, istart, istop, ii, stagout, stag1)
        ! Subroutine arguments
        INTEGER(intk), INTENT(out) :: i1
        INTEGER(intk), INTENT(inout) :: istart, istop
        INTEGER(intk), INTENT(in) :: ii, stagout, stag1

        ! Local variables
        ! none...

        IF (stagout == 0) THEN
            IF (stag1 == 0) THEN
                i1 = 0
            ELSE
                i1 = -1
            END IF
        ELSE
            IF (stag1 == 0) THEN
                i1 = 1
            ELSE
                i1 = 0
            END IF
        END IF

        IF (i1 == -1) istart = 2
        IF (i1 == 1) istop = ii - 1
    END SUBROUTINE get_stag_shift


    PURE SUBROUTINE get_idx(ip, coord, x)
        ! Subroutine arguments
        INTEGER(intk), INTENT(out) :: ip
        REAL(realk), INTENT(in) :: coord
        REAL(realk), INTENT(in), CONTIGUOUS :: x(:)

        ! Local variables
        INTEGER(intk) :: ii, i

        ii = SIZE(x)

        ip = 0
        DO i = 3, ii-1
            IF (x(i) >= coord) THEN
                ip = i-1
                EXIT
            END IF
        END DO

        ! solintxyzgrad0 access indices from ip-1 to ip+2
        IF (ip < 2 .OR. ip > ii-2) THEN
            ERROR STOP
        END IF
    END SUBROUTINE get_idx
END MODULE utils_mod
