 SUBROUTINE CREATE_IBM_DATA_STRUC

 USE GLOBAL
 USE MESH2D
 USE IBM_PARAMS
 IMPLICIT NONE

 INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IS2C_TEMP, IGC2C_TEMP
 INTEGER :: icG, jcG, iGC
 INTEGER :: iSolid

 DOUBLE PRECISION :: XmaxCyl, XminCyl, YmaxCyl, YminCyl
 DOUBLE PRECISION :: xx, yy, rad

 CALL READ_CYLINDER_INPUT

 WRITE(*,'(A)',ADVANCE='NO') ' CREATING IBM DATA STRUCTURES ...'

 !Create IBlank
 NSolid = 0; IBLANK = FLUID
 ALLOCATE(IS2C_TEMP(nxcg*nycg,2))

 IF (CYLINDER_TYPE .EQ. SQUARE) THEN

	 XmaxCyl = xcent + 0.5D0*D_cyl
	 XminCyl = xcent - 0.5D0*D_cyl
	 YmaxCyl = ycent + 0.5D0*D_cyl
	 YminCyl = ycent - 0.5D0*D_cyl

	 DO jcG = 1, nycg
		DO icG = 1, nxcg
			xx = xcg(icG,jcG); yy = ycg(icG,jcG)
			IF (xx .LE. XmaxCyl .AND. xx .GE. XminCyl .AND. &
				yy .LE. YmaxCyl .AND. yy .GE. YminCyl) THEN

				NSolid = NSolid + 1
				IS2C_TEMP(NSolid,1) = icG; IS2C_TEMP(NSolid,2) = jcG
				IBLANK(icG,jcG) = SOLID
			END IF
		END DO
	 END DO

 ELSEIF (CYLINDER_TYPE .EQ. CIRCULAR) THEN

	 DO jcG = 1, nycg
		DO icG = 1, nxcg
			xx = xcg(icG,jcG); yy = ycg(icG,jcG)
			rad = SQRT((xx-xcent)**2+(yy-ycent)**2)
			IF (rad .LE. 0.5*D_cyl) THEN
				NSolid = NSolid + 1
				IS2C_TEMP(NSolid,1) = icG; IS2C_TEMP(NSolid,2) = jcG
				IBLANK(icG,jcG) = SOLID
			END IF
		END DO
	 END DO

 END IF

 GMAP = 0

 !Copy from temporary variable
 ALLOCATE(IS2C(NSolid,2))
 IS2C(:NSolid,1) = IS2C_TEMP(:NSolid,1)
 IS2C(:NSolid,2) = IS2C_TEMP(:NSolid,2)

 !Create Ghost Cell Data structures
 NGC = 0; IGC2C_TEMP = 0
 ALLOCATE(IGC2C_TEMP(NSolid,2))

 DO iSolid = 1, NSolid
	icG = IS2C(iSolid,1); jcG = IS2C(iSolid,2)

	IF (IBLANK(icG-1,jcG) .EQ. FLUID .OR. IBLANK(icG+1,jcG) .EQ. FLUID .OR. &
		IBLANK(icG,jcG-1) .EQ. FLUID .OR. IBLANK(icG,jcG+1) .EQ. FLUID) THEN

		NGC = NGC + 1
		IGC2C_TEMP(NGC,1) = icG; IGC2C_TEMP(NGC,2) = jcG
		GMAP(icG,jcG) = 1
	END IF
 END DO

 !Copy from temporary variable
 ALLOCATE(IGC2C(NGC,2))
 IGC2C(:NGC,1) = IGC2C_TEMP(:NGC,1)
 IGC2C(:NGC,2) = IGC2C_TEMP(:NGC,2)

 !Find fluid neighbors for Ghost cells
 ALLOCATE(IGC2FNB(NGC,2))

 DO iGC = 1, NGC
	icG = IGC2C(IGC,1); jcG = IGC2C(IGC,2)

	IF 		(IBLANK(icG-1,jcG) .EQ. FLUID) THEN
		IGC2FNB(IGC,1) = icG-1; IGC2FNB(IGC,2) = jcG
	ELSEIF	(IBLANK(icG+1,jcG) .EQ. FLUID) THEN
		IGC2FNB(IGC,1) = icG+1; IGC2FNB(IGC,2) = jcG
	ELSEIF	(IBLANK(icG,jcG+1) .EQ. FLUID) THEN
		IGC2FNB(IGC,1) = icG; IGC2FNB(IGC,2) = jcG+1
	ELSEIF	(IBLANK(icG,jcG-1) .EQ. FLUID) THEN
		IGC2FNB(IGC,1) = icG; IGC2FNB(IGC,2) = jcG-1
	END IF
 END DO

 WRITE(*,'(A)') ' DONE!'
 PRINT*,''

 CALL DISP_IBM_STATS

 DEALLOCATE(IS2C_TEMP,IGC2C_TEMP)
 
 END SUBROUTINE CREATE_IBM_DATA_STRUC
 !----------------------------------------------------------------------------!


 SUBROUTINE READ_CYLINDER_INPUT

 USE GLOBAL
 USE MESH2D
 USE IBM_PARAMS
 IMPLICIT NONE

 INTEGER, PARAMETER :: fcyl = 701

 WRITE(*,'(A)',ADVANCE='NO') ' READING CYLINDER INPUT FILE ...'

 OPEN(fcyl,FILE='Cylinder_input.dat',ACTION='READ')

 READ(fcyl,*)
 READ(fcyl,*)
 READ(fcyl,*) CYLINDER_TYPE
 READ(fcyl,*)
 READ(fcyl,*)
 READ(fcyl,*) xcent, ycent, D_cyl

 CLOSE(fcyl)

 WRITE(*,'(A)') ' DONE!'
 PRINT*,''
 
 END SUBROUTINE READ_CYLINDER_INPUT
 !----------------------------------------------------------------------------!


 SUBROUTINE DISP_IBM_STATS

 USE IBM_PARAMS
 IMPLICIT NONE

 WRITE(*,*) '------------------------------------------------------------------------'
 WRITE(*,*) 'IBM STATISTICS:'
 WRITE(*,*) ''

 IF (CYLINDER_TYPE .EQ. SQUARE) THEN
	WRITE(*,'(A)') ' >> CYLINDER TYPE: SQUARE'
 ELSEIF (CYLINDER_TYPE .EQ. CIRCULAR) THEN
	WRITE(*,'(A)') ' >> CYLINDER TYPE: CIRCULAR'
 ELSE
	WRITE(*,'(A)') ' >> ERROR: UNKNOWN CYLINDER TYPE. ABORTING!'
	STOP
 END IF

 WRITE(*,'(A,I5)') ' >> NUMBER OF SOLID CELLS DETECTED: ', NSolid
 WRITE(*,'(A,I5)') ' >> NUMBER OF GHOST CELLS DETECTED: ', NGC
 WRITE(*,*) '------------------------------------------------------------------------'
 WRITE(*,*) ''

 END SUBROUTINE DISP_IBM_STATS
 !----------------------------------------------------------------------------!
