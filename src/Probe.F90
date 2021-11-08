 SUBROUTINE OPEN_PROBE_FILES

 USE GLOBAL
 IMPLICIT NONE

 INTEGER :: iProbe

 CHARACTER (LEN=50) :: fname

 WRITE(*,'(A)',ADVANCE='NO') ' OPENING PROBE FILES ...'

 DO iProbe = 1, nProbe
	WRITE(fname,'(A,I2.2,A)') 'data_probe.',iProbe,'.dat'
	OPEN(fMonBase+iProbe,FILE=TRIM(fname))
	
	WRITE(fMonBase+iProbe,*) 'VARIABLES = iter, time, u, v, p'
 END DO

 WRITE(*,'(A)') ' DONE!'
 PRINT*,''

 END SUBROUTINE OPEN_PROBE_FILES
 !----------------------------------------------------------------------------!


 SUBROUTINE WRITE_PROBE_DATA

 USE GLOBAL
 USE MESH2D
 USE FLOW_PARAMS
 IMPLICIT NONE

 INTEGER :: iProbe
 INTEGER :: icG, jcG

 DO iProbe = 1, nProbe
	icG = LOC_PROBE(iProbe,1);	jcG = LOC_PROBE(iProbe,2)

	WRITE(fMonBase+iProbe,71) iter, ctime, uc(icG,jcG), vc(icG,jcG), pc(icG,jcG)
 END DO

71 FORMAT (I6,4(2X,E16.7))

 END SUBROUTINE WRITE_PROBE_DATA
 !----------------------------------------------------------------------------!


 SUBROUTINE CLOSE_PROBE_FILES

 USE GLOBAL
 IMPLICIT NONE

 INTEGER :: iProbe

 WRITE(*,'(A)',ADVANCE='NO') ' CLOSING PROBE FILES ...'

 DO iProbe = 1, nProbe
	CLOSE(fMonBase+iProbe)
 END DO

 WRITE(*,'(A)') ' DONE!'
 PRINT*,''

 END SUBROUTINE CLOSE_PROBE_FILES
