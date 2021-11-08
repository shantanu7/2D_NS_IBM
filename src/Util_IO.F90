
 SUBROUTINE READ_INPUT

 USE GLOBAL
 USE MESH2D
 USE FLOW_PARAMS
 USE POISSON_PARAMS
 IMPLICIT NONE

 INTEGER, PARAMETER :: finp = 701
 INTEGER :: iBC, j

 WRITE(*,'(A)',ADVANCE='NO') ' READING INPUT FILE ...'

 OPEN(finp, FILE='input.inp', ACTION='READ')

 READ(finp,*)
 READ(finp,*)

 READ(finp,*)
 READ(finp,*) GRID_TYPE
 READ(finp,*)
 READ(finp,*) Xmin, Xmax, Ymin, Ymax
 READ(finp,*)
 READ(finp,*) nx, ny

 !Define additional grid count parameters
 nxc  = nx-1;  nyc = ny-1			!Number of internal cells along i-, j-
 nxg  = nx+2;  nyg = ny+2
 nxcg = nxg-1; nycg = nyg-1
 NCELL	= nxc*nyc				!Number of cells -- including Ghost cells
 NCELLG = nxcg*nycg

 READ(finp,*)
 READ(finp,*) ITERMAX, dt
 dtinv = 1.D0/dt

 READ(finp,*)
 READ(finp,*)
 READ(finp,*) Re, iIntlBody

 Reinv = 1.D0/Re

 READ(finp,*)
 READ(finp,*)
 READ(finp,*) uinit, vinit, pinit

 !Read velocity Boundary Conditions
 READ(finp,*)

 DO iBC = 1, 4
	READ(finp,*)
	READ(finp,*) BCType(iBC), BCVal(iBC)
 END DO

 !Read Pressure Boundary Conditions
 READ(finp,*)
 READ(finp,*)
 READ(finp,*) P_Solver
 READ(finp,*)
 READ(finp,*) PITERMAX, omega, P_tol
 READ(finp,*)

 DO iBC = 1, 4
	READ(finp,*)
	READ(finp,*) PBCType(iBC), PBCVal(iBC)
 END DO

 READ(finp,*)
 READ(finp,*)
 READ(finp,*) TECFREQ

 CLOSE(finp)

 WRITE(*,'(A)') ' DONE!'
 PRINT*,''

 END SUBROUTINE READ_INPUT
 !----------------------------------------------------------------------------!


 SUBROUTINE SETUP_PROBE

 USE GLOBAL
 USE MESH2D
 IMPLICIT NONE

 INTEGER, PARAMETER :: fin = 704
 INTEGER :: iProbe
 INTEGER :: icG, jcG

 WRITE(*,'(A)',ADVANCE='NO') ' READING PROBE INPUT FILE ...'

 OPEN(fin,FILE='probe_in.dat')

 READ(fin,*)
 READ(fin,*)
 READ(fin,*) nProbe

 ALLOCATE(xProbe(nProbe),yProbe(nProbe))
 ALLOCATE(uProbe(nProbe),vProbe(nProbe),pProbe(nProbe))
 ALLOCATE(LOC_PROBE(nProbe,2))

 READ(fin,*)
 READ(fin,*)
 DO iProbe = 1, nProbe
	 READ(fin,*) xProbe(iProbe), yProbe(iProbe)	
 END DO

 CLOSE(fin)

 WRITE(*,'(A)') ' DONE!'
 PRINT*,''

 WRITE(*,'(A)',ADVANCE='NO') ' FINDING PROBE LOCATIONS ...'

 DO iProbe = 1, nProbe

	DO icG = 1, nxcg
		IF ( xProbe(iProbe) .GT. xcg(icG,1) .AND. xProbe(iProbe) .LE. xcg(icg+1,1)) THEN
			LOC_PROBE(iProbe,1) = icG;
			EXIT
		END IF
	END DO
	
	DO jcG = 1, nycg
		IF ( yProbe(iProbe) .GT. ycg(1,jcg) .AND. yProbe(iProbe) .LE. ycg(1,jcg+1)) THEN
			LOC_PROBE(iProbe,2) = jcG;
			EXIT
		END IF
	END DO

 END DO

 WRITE(*,'(A)') ' DONE!'
 PRINT*,''

 DO iProbe = 1, nProbe
	WRITE(*,'(A,I1,A,2F)') ' CLOSEST CELL FOR iProbe: ', iProbe, ': ', &
	xcg(LOC_PROBE(iProbe,1),LOC_PROBE(iProbe,2)), ycg(LOC_PROBE(iProbe,1),LOC_PROBE(iProbe,2))
 END DO
 PRINT*,''

 CALL OPEN_PROBE_FILES

 END SUBROUTINE SETUP_PROBE
 !----------------------------------------------------------------------------!


 SUBROUTINE WRITE_GRID_FILE

 USE GLOBAL
 USE MESH2D
 IMPLICIT NONE

 DOUBLE PRECISION :: xUf(nx,nyc), yUf(nx,nyc)
 DOUBLE PRECISION :: xVf(nxc,ny), yVf(nxc,ny)

 INTEGER, PARAMETER :: fgrid = 801
 INTEGER :: i, j

 WRITE(*,'(A)',ADVANCE='NO') ' WRITING GRID FILE ...'

 !Create Uf grid
 DO j = 1, nyc
	xUf(:nx,j) = xv(:nx,1)
	yUf(:nx,j) = yc(1,j)
 END DO

 !Create Vf grid
 DO i = 1, nxc
	xVf(i,:ny) = xc(i,1)
	yVf(i,:ny) = yv(1,:ny)
 END DO

 OPEN(fgrid,FILE='grid_layout.dat')

 WRITE(fgrid,'(A)') ' TITLE = "STAGGERED GRID ARRANGEMENT"'
 WRITE(fgrid,'(A)') ' VARIABLES = "X", "Y"'

 !Create internal grid
 WRITE(fgrid,*) ' ZONE T = "MESH", I = ',nx,', J = ',ny,', F = POINT'
 DO j = 1, ny
	DO i = 1, nx
		WRITE(fgrid,*) xv(i,j), yv(i,j)
	END DO
 END DO

 !Create global grid
 WRITE(fgrid,*) ' ZONE T = "MESHG", I = ',nxg,', J = ',nyg,', F = POINT'
 DO j = 1, nyg
	DO i = 1, nxg
		WRITE(fgrid,*) xvg(i,j), yvg(i,j)
	END DO
 END DO

 !Write internal cell-centers
 WRITE(fgrid,*) ' ZONE T = "CC", I = ',nxc,', J = ',nyc,', F = POINT'
 DO j = 1, nyc
	DO i = 1, nxc
		WRITE(fgrid,*) xc(i,j), yc(i,j)
	END DO
 END DO

 !Write lobal cell-centers
 WRITE(fgrid,*) ' ZONE T = "CCG", I = ',nxcg,', J = ',nycg,', F = POINT'
 DO j = 1, nycg
	DO i = 1, nxcg
		WRITE(fgrid,*) xcg(i,j), ycg(i,j)
	END DO
 END DO

 !Write Face x-velocity points
 WRITE(fgrid,*) ' ZONE T = "UFC", I = ',nx,', J = ',nyc,', F = POINT'
 DO j = 1, nyc
	DO i = 1, nx
		WRITE(fgrid,*) xUf(i,j), yUf(i,j)
	END DO
 END DO

 !Write Face y-velocity points
 WRITE(fgrid,*) ' ZONE T = "VFC", I = ',nxc,', J = ',ny,', F = POINT'
 DO j = 1, ny
	DO i = 1, nxc
		WRITE(fgrid,*) xVf(i,j), yVf(i,j)
	END DO
 END DO

 CLOSE(fgrid)


 PRINT*,'DONE!'
 PRINT*,''

 END SUBROUTINE WRITE_GRID_FILE
 !----------------------------------------------------------------------------!


 SUBROUTINE WRITE_TEC_OUT_ASCII

 USE GLOBAL
 USE MESH2D
 USE FLOW_PARAMS
 USE IBM_PARAMS
 IMPLICIT NONE

 DOUBLE PRECISION :: wz(nxcg,nycg), plotQ(nxc,nyc,5)

 INTEGER :: icG, jcG, ic, jc
 INTEGER, PARAMETER :: fout = 801

 CHARACTER (LEN=50)  :: fname

 DO jcG = 2, nycg-1
	DO icG = 2, nxcg-1
		ic = icG - 1; jc = jcG - 1
		wz(icG,jcG) = (uc(icG,jcG+1)-uc(icG,jcG-1))*dycinv_NC(ic,jc) - &
					  (vc(icG+1,jcG)-vc(icG-1,jcG))*dxcinv_NC(ic,jc)
	END DO
 END DO

 DO jc = 1, nyc
	DO ic = 1, nxc
		icG = ic + 1; jcG = jc + 1
		plotQ(ic,jc,1) = (1.-IBLANK(icG,jcG))*uc(icG,jcG)
		plotQ(ic,jc,2) = (1.-IBLANK(icG,jcG))*vc(icG,jcG)
		plotQ(ic,jc,3) = (1.-IBLANK(icG,jcG))*pc(icG,jcG)
		plotQ(ic,jc,4) = (1.-IBLANK(icG,jcG))*wz(icG,jcG)
		plotQ(ic,jc,5) = IBLANK(icG,jcG)
	END DO
 END DO

 WRITE(*,'(A)',ADVANCE='NO') ' WRITING ASCII TECPLOT FILE ...'

 WRITE(fname,'(A,I5.5,A)') 'NS_out.', ITER, '.dat'
 OPEN(fout,FILE=TRIM(fname))
 WRITE(fout,*) 'VARIABLES = "X","Y","U","V","P","WZ","BL"'
 WRITE(fout,*) 'ZONE I=', nx, 'J=', ny, 'ZONETYPE=ORDERED'
 WRITE(fout,*) 'DATAPACKING=BLOCK VARLOCATION=([3,4,5,6,7]=CELLCENTERED)'

 WRITE(fout,*) xv
 WRITE(fout,*) yv
 WRITE(fout,*) plotQ(:,:,1)
 WRITE(fout,*) plotQ(:,:,2)
 WRITE(fout,*) plotQ(:,:,3)
 WRITE(fout,*) plotQ(:,:,4)
 WRITE(fout,*) plotQ(:,:,5)

 CLOSE(fout)

 PRINT*,'DONE!'

 END SUBROUTINE WRITE_TEC_OUT_ASCII
 !----------------------------------------------------------------------------!


 SUBROUTINE WRITE_TEC_OUT_BINARY

 USE GLOBAL
 USE MESH2D
 USE FLOW_PARAMS
 USE POISSON_PARAMS
 IMPLICIT NONE

 INTEGER, PARAMETER :: ftec = 801
 INTEGER :: i, j, k, iG, jG

 !Tecplot binary format related
 !-----------------------------
 INTEGER, PARAMETER :: NVARS = 5, MAX_STRING_LEN = 64 
 INTEGER, ALLOCATABLE, DIMENSION(:) :: itemp

 REAL(KIND=4), PARAMETER :: ZONEMARKER = 299.0, EOHMARKER = 357.0 
 REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: minmax
 REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: plotQ
 REAL(KIND=8) :: minU, maxU, maxV, minV, maxP, minP

 CHARACTER (LEN=100) :: string
 CHARACTER (LEN=50)  :: fname

 !Allocate memory associated with header
 ALLOCATE(minmax(NVARS,2), itemp(MAX_STRING_LEN), plotQ(nx,ny,NVARS-NDIM))
 plotQ = 0.;
 minU = 1.E10; maxU = -1.E10;
 minV = 1.E10; maxV = -1.E10;
 minP = 1.E10; maxP = -1.E10;

 DO j = 1, nyc
	DO i = 1, nxc
		iG = i + 1; jG = j + 1
		maxU = MAX(uc(iG,jG),maxU); minU = MIN(uc(iG,jG),minU)
		maxV = MAX(vc(iG,jG),maxV);	minV = MIN(vc(iG,jG),minV)
		maxP = MAX(pc(iG,jG),maxP);	minP = MIN(pc(iG,jG),minP)
	END DO
 END DO

 DO j = 1, nyc
	DO i = 1, nxc
		plotQ(i,j,1) = uc(i+1,j+1)
		plotQ(i,j,2) = vc(i+1,j+1)
		plotQ(i,j,3) = pc(i+1,j+1)
	END DO
 END DO

 WRITE(*,'(A)',ADVANCE='NO') '    WRITING BINARY TECPLOT FILE ...'

 WRITE(fname,'(A,I5.5,A)') 'NS_out.', ITER, '.plt'
 OPEN(ftec, FILE=TRIM(fname), ACCESS='STREAM', FORM='UNFORMATTED', ACTION='WRITE')

 WRITE(ftec) '#!TDV112'	!Magic number (Tecplot Version)

 !Title of file
 itemp(1) = 1			!Integer value of 1
 itemp(2) = 0			!File type (0=full, 1=grid, 2=solution)

 WRITE(ftec) itemp(:2)

 !Start writing file header
 string = '2D_Poisson_Solver_result' 	!Title of data set
 CALL WRITE_TEC_STRING(ftec, string)

 !NVARS and variable names
 WRITE(ftec) NVARS

 string = "X"
 CALL WRITE_TEC_STRING(ftec, string)

 string = "Y"
 CALL WRITE_TEC_STRING(ftec, string)

 string = "U"
 CALL WRITE_TEC_STRING(ftec, string)

 string = "V"
 CALL WRITE_TEC_STRING(ftec, string)

 string = "P"
 CALL WRITE_TEC_STRING(ftec, string)

 WRITE(ftec) ZONEMARKER

 string = 'ZONE 001'
 CALL WRITE_TEC_STRING(ftec, string)

 itemp(1) =-1		!Parent Zone
 itemp(2) = 0		!Strand ID
 WRITE(ftec) itemp(:2)

 WRITE(ftec) REAL(ctime,KIND=8) !Solution time

 !Dataset format information
 itemp(1)  = -1		!Zone color
 itemp(2)  =  0		!Zone type (0=Ordered data)
! itemp(3)  =  0		!Data packing (0=block, 1=point)
 itemp(3)  =  1		!Var location (0=nodal, 1=specify)

 WRITE(ftec) itemp(:3)

 itemp(1)  =  0		!0=nodal, 1=cellcentered
 itemp(2)  =  0		!0=nodal, 1=cellcentered
 itemp(3)  =  1		!0=nodal, 1=cellcentered
 itemp(4)  =  1		!0=nodal, 1=cellcentered
 itemp(5)  =  1		!0=nodal, 1=cellcentered
 itemp(6)  =  0		!Raw local data supplied? 0: No
 itemp(7)  =  0		!Number of misc. user-defined face conn
 itemp(8)  =  nx	!Icell dim (for i-, j-, k- type data)
 itemp(9)  =  ny	!Jcell dim (for i-, j-, k- type data)
 itemp(10) =  1		!Kcell dim (for i-, j-, k- type data)
 itemp(11) =  0		!Other auxiliary data

 WRITE(ftec) itemp(:11)

 WRITE(ftec) EOHMARKER

 !Write Zones
 WRITE(ftec) ZONEMARKER

 itemp(1:NVARS) = 2	!Variable format: 1=float, 2=double...
 itemp(NVARS+1) = 0	!Has passive variables?
 itemp(NVARS+2) = 0	!Has variable sharing?
 itemp(NVARS+3) =-1	!Zone to share connectivity (-1=no)

 WRITE(ftec) itemp(:NVARS+3)

 minmax(1,1) = Xmin; 		minmax(1,2) = Xmax
 minmax(2,1) = Ymin; 		minmax(2,2) = Ymax
 minmax(3,1) = minU; 		minmax(3,2) = maxU;
 minmax(4,1) = minV; 		minmax(4,2) = maxV;
 minmax(5,1) = minP; 		minmax(5,2) = maxP;

 DO k = 1, NVARS
 	WRITE(ftec) REAL(minmax(k,1),KIND=8)
 	WRITE(ftec) REAL(minmax(k,2),KIND=8)
 END DO

 !Start writing variables
 DO j = 1, ny
 	WRITE(ftec) REAL(xv(:nx,j),KIND=8)
 END DO

 DO j = 1, ny
 	WRITE(ftec) REAL(yv(:nx,j),KIND=8)
 END DO

 DO k = 1, NVARS-NDIM
	DO j = 1, ny
	 	WRITE(ftec) REAL(plotQ(:nx,j,k),KIND=8)
	END DO
 END DO

 CLOSE(ftec)

 PRINT*,'DONE!'
 PRINT*,''

 !Free used memory
 DEALLOCATE(minmax, itemp, plotQ) 

!-------------------------------------------------------!
 CONTAINS

 SUBROUTINE WRITE_TEC_STRING(ftec, string)

 IMPLICIT NONE

 INTEGER, INTENT(IN) :: ftec
 INTEGER, ALLOCATABLE, DIMENSION (:) :: int_str
 INTEGER :: i, strlen

 CHARACTER(LEN=100), INTENT(IN) :: string

 strlen = LEN_TRIM(string)

 ALLOCATE(int_str(strlen+1))

 DO i = 1, strlen
	int_str(i) = ICHAR(string(i:i))
 END DO
 int_str(strlen+1) = 0

 WRITE(ftec) int_str(1:strlen+1)

 DEALLOCATE(int_str)
 
 END SUBROUTINE WRITE_TEC_STRING 
 
 END SUBROUTINE WRITE_TEC_OUT_BINARY
 !----------------------------------------------------------------------------!
