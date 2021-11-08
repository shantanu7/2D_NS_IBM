 SUBROUTINE INITSETUP

 USE GLOBAL
 USE MESH2D
 USE FLOW_PARAMS
 USE POISSON_PARAMS
 IMPLICIT NONE

 INTEGER :: i, j, ic, jc


 WRITE(*,'(A)',ADVANCE='NO') ' SETTING INITIAL CONDITIONS ...'

 DO j = 2, nycg-1
	DO i = 1, nxcg
		uc(i,j) = uinit*(1.D0 - 4.D0*(ycg(i,j)/Ly)**2)
	END DO
 END DO

 DO jc = 1, nyc
	DO i = 1, nx 
		Uf(i,jc) = uinit*(1.D0 - 4.D0*(ycg(i,jc+1)/Ly)**2)
	END DO
 END DO

 vc = vinit; Vf = vinit

 usr = uc; vsr = vc

 !Pressure 
! DO j = 1, nycg
!	 DO i = 1, nxcg
!		pc(i,j) = PBCVal(2) - (xcg(i,j)-Xmin)*Reinv*(8.D0*uinit/Ly)
!	 END DO
! END DO

 pc = 0.D0

 ctime = 0.D0

 WRITE(*,'(A)') ' DONE!'
 WRITE(*,'(A)') ''

! CALL INIT_POISSON

 CALL WRITE_TEC_OUT_ASCII

 CALL SETUP_PROBE

 END SUBROUTINE INITSETUP
 !----------------------------------------------------------------------------!
