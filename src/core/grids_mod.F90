MODULE grids_mod
    USE precision_mod, ONLY: intk, realk, mglet_filename_max
    USE err_mod, ONLY: errr
    USE gridio_mod, ONLY: gridinfo_t, bcond_t, maxboconds

    IMPLICIT NONE (type, external)
    PRIVATE

    TYPE(gridinfo_t), ALLOCATABLE, TARGET :: gridinfo(:)
    TYPE(bcond_t), ALLOCATABLE, TARGET :: front(:)
    TYPE(bcond_t), ALLOCATABLE, TARGET :: back(:)
    TYPE(bcond_t), ALLOCATABLE, TARGET :: right(:)
    TYPE(bcond_t), ALLOCATABLE, TARGET :: left(:)
    TYPE(bcond_t), ALLOCATABLE, TARGET :: bottom(:)
    TYPE(bcond_t), ALLOCATABLE, TARGET :: top(:)

    REAL(realk), ALLOCATABLE :: realprms(:)
    INTEGER(intk), ALLOCATABLE :: intprms(:)

    ! Elements from former mgpar.h, colevel.h cobound.h and setmpi_mod
    INTEGER(intk), PROTECTED :: ngrid
    INTEGER(intk), PROTECTED :: minlevel
    INTEGER(intk), PROTECTED :: maxlevel
    INTEGER(intk), PROTECTED :: maxgrdsoflvl

    INTEGER(intk), ALLOCATABLE, PROTECTED :: noflevel(:), igrdoflevel(:, :)

    INTEGER(intk), PROTECTED :: nmygrids
    INTEGER(intk), ALLOCATABLE, PROTECTED :: mygrids(:)
    INTEGER(intk), ALLOCATABLE, PROTECTED :: nmygridslvl(:)
    INTEGER(intk), ALLOCATABLE, PROTECTED :: mygridslvl(:, :)

    ! From cobound.h
    INTEGER(intk), ALLOCATABLE, PROTECTED :: nboconds(:, :)
    INTEGER(intk), ALLOCATABLE, PROTECTED :: itypboconds(:, :, :)

    ! From compi.h
    INTEGER(intk), ALLOCATABLE, PROTECTED :: idprocofgrd(:)

    INTERFACE get_mgbasb
        MODULE PROCEDURE :: get_mgbasb1, get_mgbasb2
    END INTERFACE get_mgbasb

    ! Public subroutines
    PUBLIC :: init_grids, finish_grids, get_bbox, get_gradpxflag, get_bcprms, &
        get_imygrid, get_mgdims, get_level, iposition, jposition, kposition, &
        iparent, get_neighbours, level, get_kk, get_jj, get_ii, get_mgbasb, &
        get_bc_ctyp, get_gridvolume, get_nbcprms

    ! Public data arrays
    PUBLIC :: ngrid, minlevel, maxlevel, maxgrdsoflvl, noflevel, igrdoflevel, &
        nmygrids, mygrids, nmygridslvl, mygridslvl, nboconds, itypboconds, &
        idprocofgrd

CONTAINS

    SUBROUTINE init_grids()
        USE HDF5
        USE gridio_mod, ONLY: read_gridinfo, read_bcondinfo
        USE hdf5common_mod, ONLY: hdf5common_open, hdf5common_close
        USE comms_mod, ONLY: numprocs
        USE fort7_mod, ONLY: fort7

        ! Local variables
        INTEGER(HID_T) :: file_id
        CHARACTER(len=mglet_filename_max) :: filename

        CALL fort7%get_value("/io/grids", filename, "grids.h5")
        CALL hdf5common_open(filename, "r", file_id)

        CALL read_gridinfo(file_id, gridinfo, realprms, intprms, ngrid)

        ALLOCATE(front(ngrid))
        CALL read_bcondinfo(file_id, "FRONT", front)

        ALLOCATE(back(ngrid))
        CALL read_bcondinfo(file_id, "BACK", back)

        ALLOCATE(right(ngrid))
        CALL read_bcondinfo(file_id, "RIGHT", right)

        ALLOCATE(left(ngrid))
        CALL read_bcondinfo(file_id, "LEFT", left)

        ALLOCATE(bottom(ngrid))
        CALL read_bcondinfo(file_id, "BOTTOM", bottom)

        ALLOCATE(top(ngrid))
        CALL read_bcondinfo(file_id, "TOP", top)

        CALL hdf5common_close(file_id)

        ! Sanity check
        IF (ngrid < numprocs) THEN
            WRITE(*,*) "ngrid:", ngrid
            WRITE(*,*) "Less than numprocs:", numprocs
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Set colevel.h, minlevel, maxlevel, noflevel, igrdoflevel
        CALL setcolevel()

        ! Distribute grids between processes, idprocofgrd
        CALL setmpi()

        ! Set boundary conditions, pointers
        CALL init_gridstructure()
    END SUBROUTINE init_grids


    SUBROUTINE finish_grids()
        DEALLOCATE(gridinfo)
        DEALLOCATE(front)
        DEALLOCATE(back)
        DEALLOCATE(right)
        DEALLOCATE(left)
        DEALLOCATE(bottom)
        DEALLOCATE(top)
        DEALLOCATE(realprms)
        DEALLOCATE(intprms)

        DEALLOCATE(noflevel)
        DEALLOCATE(igrdoflevel)
        DEALLOCATE(mygrids)
        DEALLOCATE(nmygridslvl)
        DEALLOCATE(mygridslvl)
        DEALLOCATE(nboconds)
        DEALLOCATE(itypboconds)
        DEALLOCATE(idprocofgrd)

        ngrid = 0
        nmygrids = 0
        minlevel = 0
        maxlevel = 0
        maxgrdsoflvl = 0
    END SUBROUTINE finish_grids


    SUBROUTINE setcolevel()
        ! Subroutine arguments
        ! none...

        ! Local variables
        INTEGER(intk) :: igrid, level

        ! Minimum and maximum levels
        minlevel = gridinfo(1)%level
        maxlevel = gridinfo(ngrid)%level
        ALLOCATE(noflevel(minlevel:maxlevel))
        noflevel = 0

        ! Determine number of grids per level, maximum number of grids per
        ! level
        maxgrdsoflvl = 0
        DO igrid = 1, ngrid
            level = gridinfo(igrid)%level
            noflevel(level) = noflevel(level) + 1
            maxgrdsoflvl = MAX(maxgrdsoflvl, noflevel(level))
        END DO

        ! Allocate and set igrdoflevel, re-set noflevel (needed as counter)
        ALLOCATE(igrdoflevel(maxgrdsoflvl, minlevel:maxlevel))
        noflevel = 0
        igrdoflevel = 0

        DO igrid = 1, ngrid
            level = gridinfo(igrid)%level
            noflevel(level) = noflevel(level) + 1
            igrdoflevel(noflevel(level), level) = igrid
        END DO
    END SUBROUTINE setcolevel


    SUBROUTINE init_gridstructure()
        !USE pointer_mod, ONLY: init_pointers, set_pointer

        INTEGER(intk) :: igrid

        ALLOCATE(nboconds(6, ngrid))
        ALLOCATE(itypboconds(maxboconds, 6, ngrid))

        nboconds = 0.0
        itypboconds = 0

        !CALL init_pointers(ngrid)
        DO igrid = 1, ngrid
            ! Set boundary conditions in cobound.h
            CALL setcobound(igrid, 1, front)
            CALL setcobound(igrid, 2, back)
            CALL setcobound(igrid, 3, right)
            CALL setcobound(igrid, 4, left)
            CALL setcobound(igrid, 5, bottom)
            CALL setcobound(igrid, 6, top)

            ! Set pointers for 3D, 2D, 1D storage
            !CALL set_pointer(igrid, kk, jj, ii)
        END DO
    END SUBROUTINE init_gridstructure


    SUBROUTINE setcobound(igrid, idir, face)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(in) :: idir
        TYPE(bcond_t), INTENT(in) :: face(:)

        ! Local variables
        INTEGER(intk) :: i, nbocd
        CHARACTER(len=SIZE(face(1)%type, 1)) :: type

        nbocd = face(igrid)%nbocd
        DO i = 1, nbocd
            type = TRANSFER(face(igrid)%type(:, i), type)
            CALL setcobone(igrid, idir, i, TRIM(type))
        END DO
    END SUBROUTINE setcobound


    SUBROUTINE setcobone(igrid, idir, ibocond, ctyp)
        INTEGER, INTENT(IN) :: igrid, idir, ibocond
        CHARACTER(LEN=*), INTENT(IN) :: ctyp

        nboconds(idir, igrid) = MAX(nboconds(idir, igrid), ibocond)

        IF (ctyp == 'FIX') THEN
            itypboconds(ibocond, idir, igrid) = 2
        ELSEIF (ctyp == 'OP1') THEN
            itypboconds(ibocond, idir, igrid) = 3
        ELSEIF (ctyp == 'NOS') THEN
            itypboconds(ibocond, idir, igrid) = 5
        ELSEIF (ctyp == 'SLI') THEN
            itypboconds(ibocond, idir, igrid) = 6
        ELSEIF (ctyp == 'CON') THEN
            itypboconds(ibocond, idir, igrid) = 7
        ELSEIF (ctyp == 'PAR') THEN
            itypboconds(ibocond, idir, igrid) = 8
        ELSEIF (ctyp == 'GRA') THEN
            itypboconds(ibocond, idir, igrid) = 15
        ELSEIF (ctyp == 'SCA') THEN
            itypboconds(ibocond, idir, igrid) = 16
        ELSEIF (ctyp == 'REF') THEN
            itypboconds(ibocond, idir, igrid) = 17
        ELSEIF (ctyp == 'NRE') THEN
            itypboconds(ibocond, idir, igrid) = 18
        ELSEIF (ctyp == 'CO1') THEN
            IF (idir == 2 .OR. idir == 4 .OR. idir == 6) THEN
                WRITE(6,*) 'SETCOBONE: CO1 NOT ALLOWED IN IDIR=',IDIR
                CALL errr(__FILE__, __LINE__)
            END IF
            itypboconds(ibocond, idir, igrid) = 19
        ELSEIF (CTYP .EQ. 'NIX') THEN
            itypboconds(ibocond, idir, igrid) = 99
        ELSE
            itypboconds(ibocond, idir, igrid) = -HUGE(1_intk)
        END IF
    END SUBROUTINE setcobone


    SUBROUTINE setmpi()
        USE comms_mod, ONLY: myid

        INTEGER(intk) :: ilevel, igrid, i

        ALLOCATE(idprocofgrd(ngrid))
        CALL dist_grids()

        ! Set number of grids on current process
        ALLOCATE(nMyGridsLvl(minlevel:maxlevel))
        nMyGrids = 0
        nMyGridsLvl = 0
        DO ilevel = minlevel, maxlevel
            DO i = 1, noflevel(ilevel)
                igrid = igrdoflevel(i, ilevel)
                IF (myid == idprocofgrd(igrid)) THEN
                    nMyGrids = nMyGrids + 1
                    nMyGridsLvl(ilevel) = nMyGridsLvl(ilevel) + 1
                END IF
            END DO
        END DO

        ALLOCATE(myGrids(nMyGrids))
        ALLOCATE(myGridsLvl(MAXVAL(nMyGridsLvl), minlevel:maxlevel))
        nMyGrids = 0
        myGrids = 0
        nMyGridsLvl = 0
        myGridsLvl = 0
        DO ilevel = minlevel, maxlevel
            DO i = 1, noflevel(ilevel)
                igrid = igrdoflevel(i, ilevel)
                IF (myid == idprocofgrd(igrid)) THEN
                    nMyGrids = nMyGrids + 1
                    myGrids(nMyGrids) = igrid

                    nMyGridsLvl(ilevel) = nMyGridsLvl(ilevel) + 1
                    myGridsLvl(nMyGridsLvl(ilevel), ilevel) = igrid
                END IF
            END DO
        END DO

        IF (myid == 0) CALL setmpi_info()
    END SUBROUTINE setmpi


    SUBROUTINE dist_grids()
        USE comms_mod, ONLY: numprocs

        INTEGER(intk) :: ilevel, iproc, igrid
        INTEGER(intk) :: i, n, restProc, rest
        INTEGER(intk) :: grdProc

        INTEGER(intk), ALLOCATABLE :: nGrdsOfProc(:)

        ALLOCATE(nGrdsOfProc(0:numprocs-1))

        restProc = 0
        DO ilevel = minlevel, maxlevel
            ! Reset number of grids per process
            nGrdsOfProc = 0

            ! STEP 1: Distribute even part of grids
            n = noflevel(ilevel)/numprocs
            nGrdsOfProc = n

            ! STEP 2: Distribute rest grids
            rest = MOD(noflevel(ilevel), numprocs)
            DO i = 1, rest
                nGrdsOfProc(restProc) = nGrdsOfProc(restProc) + 1
                restProc = restProc + 1
                IF (restProc > numprocs - 1) restProc = 0
            END DO

            ! Hand out specific grids to processes
            iproc = 0     ! Rank we are handing out grids to
            grdProc = 0   ! Grids handed out to current rank
            DO i = 1, noflevel(ilevel)
                igrid = igrdoflevel(i, ilevel)

                DO WHILE (nGrdsOfProc(iproc) < 1)
                    iproc = iproc + 1
                END DO
                IF (iproc > numprocs - 1) THEN
                    CALL errr(__FILE__, __LINE__)
                END IF
                IF (grdProc < nGrdsOfProc(iproc)) THEN
                    idprocofgrd(igrid) = iproc
                    grdproc = grdproc + 1
                END IF
                IF (grdProc == nGrdsOfProc(iproc)) THEN
                    iproc = iproc + 1
                    grdproc = 0
                END IF
            END DO
        END DO
    END SUBROUTINE dist_grids


    SUBROUTINE setmpi_info()
        USE comms_mod, ONLY: numprocs

        INTEGER(intk) :: ilevel, iproc, igrid, i, grdsum
        INTEGER(intk) :: nGrdsOfProcTmp(minlevel:maxlevel)
        INTEGER(intk), ALLOCATABLE :: nGrdsOfProc(:,:)
        INTEGER(intk), ALLOCATABLE :: lvl(:)
        INTEGER(intk) :: nlvl, width

        CHARACTER(LEN=256)  :: fmt

        ALLOCATE(nGrdsOfProc(0:numprocs-1, minlevel:maxlevel))
        ALLOCATE(lvl(minlevel:maxlevel))
        nGrdsOfProc = 0
        lvl = 0

        ! Count number of grids per process and level
        DO ilevel = minlevel,maxlevel
            DO i = 1, noflevel(ilevel)
                igrid = igrdoflevel(i, ilevel)
                iproc = idprocofgrd(igrid)
                nGrdsOfProc(iproc, ilevel) = nGrdsOfProc(iproc, ilevel) + 1
            END DO
            lvl(ilevel) = ilevel
        END DO

        nlvl = maxlevel - minlevel + 1
        width = MAX(28, 13 + 7*nlvl + 7)

        WRITE(*, '("GRIDS PER LEVEL AND PROCESS:")')
        WRITE(*, '(A)') REPEAT("=", width)

        fmt = '(1X, "Level: ", 5X, '//REPEAT('I6, 1X, ', nlvl)//'3X, "Sum")'
        WRITE(*, fmt) lvl
        WRITE(*, '(A)') REPEAT("-", width)

        fmt = '(1X, "Rank ", I5, ":", 1X, '//REPEAT('I6, 1X, ', nlvl+1)//")"
        DO iproc = 0, numprocs-1
            grdsum = 0
            DO ilevel = minlevel, maxlevel
                grdsum = grdsum + nGrdsOfProc(iproc, ilevel)
                nGrdsOfProcTmp(ilevel) = nGrdsOfProc(iproc, ilevel)
            END DO
            WRITE(*,fmt) iproc, nGrdsOfProcTmp(:), grdsum
        END DO
        WRITE(*, '(A)') REPEAT("-", width)

        ! Calculate and print number of grids per level
        DO ilevel=minlevel,maxlevel
            lvl(ilevel) = SUM(nGrdsOfProc(:, ilevel))
        END DO
        fmt = '(1X, "Sum:", 8X, '//REPEAT('I6, 1X, ', nlvl+1)//")"
        WRITE(*,fmt) lvl, SUM(nGrdsOfProc(:, :))
        WRITE(*, '(A)') REPEAT("=", width)
        WRITE(*,'()')

        DEALLOCATE(nGrdsOfProc)
        DEALLOCATE(lvl)
    END SUBROUTINE setmpi_info


    SUBROUTINE get_mgdims(kk, jj, ii, igrid)
        INTEGER(intk), INTENT(OUT) :: kk, jj, ii
        INTEGER(intk), INTENT(IN) :: igrid

#ifdef _MGLET_DEBUG_
        IF (igrid > ngrid .OR. igrid < 0) THEN
            WRITE(*,*) "Invalid igrid: ", igrid
            WRITE(*,*) "ngrid: ", ngrid
            CALL errr(__FILE__, __LINE__)
        END IF
#endif

        kk = gridinfo(igrid)%kk
        jj = gridinfo(igrid)%jj
        ii = gridinfo(igrid)%ii
    END SUBROUTINE get_mgdims


    SUBROUTINE get_imygrid(imygrid, igrid)
      INTEGER(intk), INTENT(in) :: igrid
      INTEGER(intk), INTENT(out) :: imygrid

      LOGICAL :: found

      found = .FALSE.
      DO imygrid = 1, nmygrids
         IF (mygrids(imygrid) == igrid) THEN
             found = .TRUE.
             EXIT
         END IF
      END DO

      IF (.NOT. found) CALL errr(__FILE__, __LINE__)
    END SUBROUTINE get_imygrid


    SUBROUTINE get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)
        REAL(realk), INTENT(OUT) :: minx, maxx, miny, maxy, minz, maxz
        INTEGER(intk), INTENT(IN) :: igrid

#ifdef _MGLET_DEBUG_
        IF (igrid > ngrid .OR. igrid < 0) THEN
            WRITE(*,*) "Invalid igrid: ", igrid
            WRITE(*,*) "ngrid: ", ngrid
            CALL errr(__FILE__, __LINE__)
        END IF
#endif

        minx = gridinfo(igrid)%bbox(1)
        maxx = gridinfo(igrid)%bbox(2)
        miny = gridinfo(igrid)%bbox(3)
        maxy = gridinfo(igrid)%bbox(4)
        minz = gridinfo(igrid)%bbox(5)
        maxz = gridinfo(igrid)%bbox(6)
    END SUBROUTINE get_bbox


    SUBROUTINE get_gridvolume(volume, igrid)
        ! Subroutine arguments
        REAL(realk), INTENT(out) :: volume
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz

        CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)
        volume = (maxx - minx)*(maxy - miny)*(maxz - minz)
    END SUBROUTINE get_gridvolume


    SUBROUTINE get_gradpxflag(flag, igrid)
        USE simdfunctions_mod, ONLY: l_to_i

        INTEGER(intk), INTENT(OUT) :: flag
        INTEGER(intk), INTENT(IN) :: igrid

#ifdef _MGLET_DEBUG_
        IF (igrid > ngrid .OR. igrid < 0) THEN
            WRITE(*,*) "Invalid igrid: ", igrid
            WRITE(*,*) "ngrid: ", ngrid
            CALL errr(__FILE__, __LINE__)
        END IF
#endif

        flag = l_to_i(BTEST(gridinfo(igrid)%flags, 0))
    END SUBROUTINE get_gradpxflag


    SUBROUTINE get_bcprms(params, igrid, iface, ibocd)
        USE gridio_mod, ONLY: maxboconds

        ! Subroutine arguments
        CLASS(*), INTENT(INOUT) :: params(:)
        INTEGER(intk), INTENT(IN) :: igrid
        INTEGER(intk), INTENT(IN) :: iface
        INTEGER(intk), INTENT(IN) :: ibocd

        ! Local variables
        INTEGER(intk) :: offset, length, nbocd, nparams, nparams_tot
        TYPE(bcond_t), POINTER :: facearr(:)

        IF (igrid > ngrid .OR. igrid < 0) THEN
            WRITE(*,*) "Invalid igrid: ", igrid
            WRITE(*,*) "ngrid: ", ngrid
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (ibocd > maxboconds) THEN
            WRITE(*,*) "Invalid ibocd: ", ibocd
            WRITE(*,*) "maxboconds: ", maxboconds
            CALL errr(__FILE__, __LINE__)
        END IF

        SELECT CASE(iface)
        CASE (1)
            facearr => front
        CASE (2)
            facearr => back
        CASE (3)
            facearr => right
        CASE (4)
            facearr => left
        CASE (5)
            facearr => bottom
        CASE (6)
            facearr => top
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        nbocd = facearr(igrid)%nbocd
        IF (ibocd > nbocd) THEN
            WRITE(*,*) "Invalid ibocd: ", ibocd
            WRITE(*,*) "nbocd: ", nbocd
            WRITE(*,*) "igrid: ", igrid
            WRITE(*,*) "iface: ", iface
            CALL errr(__FILE__, __LINE__)
        END IF

        SELECT TYPE (params)
        TYPE IS (REAL(realk))
            offset = facearr(igrid)%realprm(1, ibocd)
            length = facearr(igrid)%realprm(2, ibocd)
            nparams = SIZE(params)
            nparams_tot = SIZE(realprms)
        TYPE IS (INTEGER(intk))
            offset = facearr(igrid)%intprm(1, ibocd)
            length = facearr(igrid)%intprm(2, ibocd)
            nparams = SIZE(params)
            nparams_tot = SIZE(intprms)
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        IF (length /= nparams) THEN
            WRITE(*,*) "Invalid length: ", length
            WRITE(*,*) "nparams: ", nparams
            WRITE(*,*) "igrid: ", igrid
            WRITE(*,*) "iface: ", iface
            CALL errr(__FILE__, __LINE__)
        END IF

        ! It is OK to ask for zero parameters, then we just return here...
        ! In this case offset == 0 is also accpeted without error...
        IF (length == 0) RETURN

        IF (offset < 1 .OR. offset+length-1 > nparams_tot) THEN
            WRITE(*,*) "Invalid offset, length: ", offset, length
            WRITE(*,*) "nparams_tot: ", nparams_tot
            WRITE(*,*) "igrid: ", igrid
            WRITE(*,*) "iface: ", iface
            CALL errr(__FILE__, __LINE__)
        END IF

        SELECT TYPE (params)
        TYPE IS (REAL(realk))
            params(1:length) = realprms(offset:offset+length-1)
        TYPE IS (INTEGER(intk))
            params(1:length) = intprms(offset:offset+length-1)
        CLASS DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE get_bcprms


    SUBROUTINE get_nbcprms(igrid, iface, ibocd, nreal, ninteger)
        ! Subroutine arguments
        INTEGER(intk), INTENT(IN) :: igrid
        INTEGER(intk), INTENT(IN) :: iface
        INTEGER(intk), INTENT(IN) :: ibocd
        INTEGER(intk), INTENT(OUT), OPTIONAL :: nreal, ninteger

        ! Local variables
        INTEGER(intk) :: nbocd
        TYPE(bcond_t), POINTER :: facearr(:)
        INTEGER(intk) :: ngrid

        ngrid = SIZE(gridinfo)
        IF (igrid > ngrid) THEN
            WRITE(*,*) "Invalid igrid: ", igrid
            WRITE(*,*) "ngrid: ", ngrid
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (ibocd > maxboconds) THEN
            WRITE(*,*) "Invalid ibocd: ", ibocd
            WRITE(*,*) "maxboconds: ", maxboconds
            CALL errr(__FILE__, __LINE__)
        END IF

        SELECT CASE(iface)
        CASE (1)
            facearr => front
        CASE (2)
            facearr => back
        CASE (3)
            facearr => right
        CASE (4)
            facearr => left
        CASE (5)
            facearr => bottom
        CASE (6)
            facearr => top
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        nbocd = facearr(igrid)%nbocd
        IF (ibocd > nbocd) THEN
            WRITE(*,*) "Invalid ibocd: ", ibocd
            WRITE(*,*) "nbocd: ", nbocd
            WRITE(*,*) "igrid: ", igrid
            WRITE(*,*) "iface: ", iface
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (PRESENT(nreal)) THEN
            nreal = facearr(igrid)%realprm(2, ibocd)
        END IF

        IF (PRESENT(ninteger)) THEN
            ninteger = facearr(igrid)%intprm(2, ibocd)
        END IF
    END SUBROUTINE get_nbcprms


    SUBROUTINE get_level(level, igrid)
        INTEGER(intk), INTENT(OUT) :: level
        INTEGER(intk), INTENT(IN) :: igrid

#ifdef _MGLET_DEBUG_
        IF (igrid > ngrid .OR. igrid < 0) THEN
            WRITE(*,*) "Invalid igrid: ", igrid
            WRITE(*,*) "ngrid: ", ngrid
            CALL errr(__FILE__, __LINE__)
        END IF
#endif

        level = gridinfo(igrid)%level
    END SUBROUTINE get_level


    INTEGER(intk) FUNCTION iposition(igrid)
        INTEGER(intk), INTENT(IN) :: igrid

#ifdef _MGLET_DEBUG_
        IF (igrid > ngrid .OR. igrid < 0) THEN
            WRITE(*,*) "Invalid igrid: ", igrid
            WRITE(*,*) "ngrid: ", ngrid
            CALL errr(__FILE__, __LINE__)
        END IF
#endif

        iposition = gridinfo(igrid)%iposition
    END FUNCTION iposition


    INTEGER(intk) FUNCTION jposition(igrid)
        INTEGER(intk), INTENT(IN) :: igrid

#ifdef _MGLET_DEBUG_
        IF (igrid > ngrid .OR. igrid < 0) THEN
            WRITE(*,*) "Invalid igrid: ", igrid
            WRITE(*,*) "ngrid: ", ngrid
            CALL errr(__FILE__, __LINE__)
        END IF
#endif

        jposition = gridinfo(igrid)%jposition
    END FUNCTION jposition


    INTEGER(intk) FUNCTION kposition(igrid)
        INTEGER(intk), INTENT(IN) :: igrid

#ifdef _MGLET_DEBUG_
        IF (igrid > ngrid .OR. igrid < 0) THEN
            WRITE(*,*) "Invalid igrid: ", igrid
            WRITE(*,*) "ngrid: ", ngrid
            CALL errr(__FILE__, __LINE__)
        END IF
#endif

        kposition = gridinfo(igrid)%kposition
    END FUNCTION kposition


    INTEGER(intk) FUNCTION iparent(igrid)
        INTEGER(intk), INTENT(IN) :: igrid

#ifdef _MGLET_DEBUG_
        IF (igrid > ngrid .OR. igrid < 0) THEN
            WRITE(*,*) "Invalid igrid: ", igrid
            WRITE(*,*) "ngrid: ", ngrid
            CALL errr(__FILE__, __LINE__)
        END IF
#endif

        iparent = gridinfo(igrid)%iparent
    END FUNCTION iparent


    INTEGER(intk) FUNCTION level(igrid)
        INTEGER(intk), INTENT(IN) :: igrid

#ifdef _MGLET_DEBUG_
        IF (igrid > ngrid .OR. igrid < 0) THEN
            WRITE(*,*) "Invalid igrid: ", igrid
            WRITE(*,*) "ngrid: ", ngrid
            CALL errr(__FILE__, __LINE__)
        END IF
#endif

        level = gridinfo(igrid)%level
    END FUNCTION level


    SUBROUTINE get_neighbours(neighbours, igrid)
        INTEGER(intk), INTENT(OUT) :: neighbours(26)
        INTEGER(intk), INTENT(IN) :: igrid

#ifdef _MGLET_DEBUG_
        IF (igrid > ngrid .OR. igrid < 0) THEN
            WRITE(*,*) "Invalid igrid: ", igrid
            WRITE(*,*) "ngrid: ", ngrid
            CALL errr(__FILE__, __LINE__)
        END IF
#endif

        neighbours = gridinfo(igrid)%nbrgrid
    END SUBROUTINE get_neighbours


    SUBROUTINE get_kk(kk, igrid)
        INTEGER(intk), INTENT(OUT) :: kk
        INTEGER(intk), INTENT(IN) :: igrid

#ifdef _MGLET_DEBUG_
        IF (igrid > ngrid .OR. igrid < 0) THEN
            WRITE(*,*) "Invalid igrid: ", igrid
            WRITE(*,*) "ngrid: ", ngrid
            CALL errr(__FILE__, __LINE__)
        END IF
#endif

        kk = gridinfo(igrid)%kk
    END SUBROUTINE get_kk


    SUBROUTINE get_jj(jj, igrid)
        INTEGER(intk), INTENT(OUT) :: jj
        INTEGER(intk), INTENT(IN) :: igrid

#ifdef _MGLET_DEBUG_
        IF (igrid > ngrid .OR. igrid < 0) THEN
            WRITE(*,*) "Invalid igrid: ", igrid
            WRITE(*,*) "ngrid: ", ngrid
            CALL errr(__FILE__, __LINE__)
        END IF
#endif

        jj = gridinfo(igrid)%jj
    END SUBROUTINE get_jj


    SUBROUTINE get_ii(ii, igrid)
        INTEGER(intk), INTENT(OUT) :: ii
        INTEGER(intk), INTENT(IN) :: igrid

#ifdef _MGLET_DEBUG_
        IF (igrid > ngrid .OR. igrid < 0) THEN
            WRITE(*,*) "Invalid igrid: ", igrid
            WRITE(*,*) "ngrid: ", ngrid
            CALL errr(__FILE__, __LINE__)
        END IF
#endif

        ii = gridinfo(igrid)%ii
    END SUBROUTINE get_ii


    SUBROUTINE get_mgbasb1(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)
        INTEGER(intk), INTENT(out) :: nfro, nbac, nrgt, nlft, nbot, ntop
        INTEGER(intk), INTENT(in) :: igrid

#ifdef _MGLET_DEBUG_
        IF (igrid > ngrid .OR. igrid < 0) THEN
            WRITE(*,*) "Invalid igrid: ", igrid
            WRITE(*,*) "ngrid: ", ngrid
            CALL errr(__FILE__, __LINE__)
        END IF
#endif

        nfro = itypboconds(1, 1, igrid)
        nbac = itypboconds(1, 2, igrid)
        nrgt = itypboconds(1, 3, igrid)
        nlft = itypboconds(1, 4, igrid)
        nbot = itypboconds(1, 5, igrid)
        ntop = itypboconds(1, 6, igrid)
    END SUBROUTINE get_mgbasb1


    SUBROUTINE get_mgbasb2(bconds, igrid)
        INTEGER(intk), INTENT(out) :: bconds(6)
        INTEGER(intk), INTENT(in) :: igrid

#ifdef _MGLET_DEBUG_
        IF (igrid > ngrid .OR. igrid < 0) THEN
            WRITE(*,*) "Invalid igrid: ", igrid
            WRITE(*,*) "ngrid: ", ngrid
            CALL errr(__FILE__, __LINE__)
        END IF
#endif

        bconds(1) = itypboconds(1, 1, igrid)
        bconds(2) = itypboconds(1, 2, igrid)
        bconds(3) = itypboconds(1, 3, igrid)
        bconds(4) = itypboconds(1, 4, igrid)
        bconds(5) = itypboconds(1, 5, igrid)
        bconds(6) = itypboconds(1, 6, igrid)
    END SUBROUTINE get_mgbasb2


    SUBROUTINE get_bc_ctyp(ctyp, ibocd, iface, igrid)
        ! Subroutine arguments
        CHARACTER(len=*), INTENT(out) :: ctyp
        INTEGER(intk), INTENT(in) :: ibocd
        INTEGER(intk), INTENT(in) :: iface
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: nchar

        ! Initialize INTENT(out)
        ctyp = ''

        IF (igrid > ngrid .OR. igrid < 0) THEN
            WRITE(*,*) "Invalid igrid: ", igrid
            WRITE(*,*) "ngrid: ", ngrid
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (iface > 6 .OR. iface < 0) THEN
            WRITE(*,*) "Invalid iface: ", iface
            CALL errr(__FILE__, __LINE__)
        END IF

        IF (ibocd > nboconds(iface, igrid) .OR. ibocd < 0) THEN
            WRITE(*,*) "Invalid ibocd: ", ibocd
            WRITE(*,*) "nboconds: ", nboconds(iface, igrid)
            CALL errr(__FILE__, __LINE__)
        END IF

        SELECT CASE (iface)
        CASE (1)
            nchar = SIZE(front(igrid)%type(:, ibocd))
            IF (nchar > LEN(ctyp)) CALL errr(__FILE__, __LINE__)
            ctyp = TRANSFER(front(igrid)%type(:, ibocd), ctyp)
        CASE (2)
            nchar = SIZE(back(igrid)%type(:, ibocd))
            IF (nchar > LEN(ctyp)) CALL errr(__FILE__, __LINE__)
            ctyp = TRANSFER(back(igrid)%type(:, ibocd), ctyp)
        CASE (3)
            nchar = SIZE(right(igrid)%type(:, ibocd))
            IF (nchar > LEN(ctyp)) CALL errr(__FILE__, __LINE__)
            ctyp = TRANSFER(right(igrid)%type(:, ibocd), ctyp)
        CASE (4)
            nchar = SIZE(left(igrid)%type(:, ibocd))
            IF (nchar > LEN(ctyp)) CALL errr(__FILE__, __LINE__)
            ctyp = TRANSFER(left(igrid)%type(:, ibocd), ctyp)
        CASE (5)
            nchar = SIZE(bottom(igrid)%type(:, ibocd))
            IF (nchar > LEN(ctyp)) CALL errr(__FILE__, __LINE__)
            ctyp = TRANSFER(bottom(igrid)%type(:, ibocd), ctyp)
        CASE (6)
            nchar = SIZE(top(igrid)%type(:, ibocd))
            IF (nchar > LEN(ctyp)) CALL errr(__FILE__, __LINE__)
            ctyp = TRANSFER(top(igrid)%type(:, ibocd), ctyp)
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT
    END SUBROUTINE get_bc_ctyp
END MODULE grids_mod
