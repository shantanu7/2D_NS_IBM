 SUBROUTINE MAKE_GRID

 USE GLOBAL
 USE MESH2D
 IMPLICIT NONE

 IF     (GRID_TYPE .EQ. UNIFORM) THEN
	CALL GRID_UNIFORM
 ELSEIF (GRID_TYPE .EQ. NON_UNIFORM) THEN
	CALL GRID_NON_UNIFORM
 END IF

 CALL DISP_GRID_STATS

 IF (iIntlBody) CALL CREATE_IBM_DATA_STRUC

 CALL CREATE_GRID_DATA_STRUC

 CALL WRITE_GRID_FILE

 END SUBROUTINE MAKE_GRID
 !----------------------------------------------------------------------------!


 SUBROUTINE GRID_UNIFORM

 USE GLOBAL
 USE MESH2D
 IMPLICIT NONE

 INTEGER :: i, j, ICELL

 WRITE(*,'(A)',ADVANCE='NO') ' MAKING UNIFORM GRID ...'

 Lx = Xmax-Xmin; Ly = Ymax-Ymin

 dx = Lx/REAL(nxc,KIND=8); dy = Ly/REAL(nyc,KIND=8)

 !Create grid points
 DO j = 1, ny
 	DO i = 1, nx
		xv(i,j) = Xmin + (i-1)*dx
		yv(i,j) = Ymin + (j-1)*dy
	END DO
 END DO

 !Create cell-centers
 DO j = 1, nyc
	DO i = 1, nxc
		xc(i,j) = 0.5D0*(xv(i,j) + xv(i+1,j))
		yc(i,j) = 0.5D0*(yv(i,j) + yv(i,j+1))
	END DO
 END DO

 !Create grid points with ghost points
 DO j = 1, nyg
 	DO i = 1, nxg
		xvg(i,j) = Xmin + (i-2)*dx
		yvg(i,j) = Ymin + (j-2)*dy
	END DO
 END DO

 !Create cell-centers with ghost cells
 DO j = 1, nycg
	DO i = 1, nxcg
		xcg(i,j) = 0.5D0*(xvg(i,j) + xvg(i+1,j))
		ycg(i,j) = 0.5D0*(yvg(i,j) + yvg(i,j+1))
	END DO
 END DO

 WRITE(*,'(A)') ' DONE!'
 WRITE(*,'(A)') ''

 END SUBROUTINE GRID_UNIFORM
 !----------------------------------------------------------------------------!


 SUBROUTINE GRID_NON_UNIFORM

 USE GLOBAL
 USE MESH2D
 IMPLICIT NONE

 INTEGER, PARAMETER :: fgrid = 701
 INTEGER :: i, j, dummy

 DOUBLE PRECISION :: xxv(nx), yyv(ny)

 WRITE(*,'(A)',ADVANCE='NO') ' READING NON-UNIFORM GRID ...'

 !Read domain vertices
 OPEN(fgrid,FILE='xgrid.dat',ACTION='READ')
 DO i = 1, nx
	READ(fgrid,*)dummy, xxv(i)
	xv(i,:) = xxv(i)
 END DO
 CLOSE(fgrid)

 OPEN(fgrid,FILE='ygrid.dat',ACTION='READ')
 DO j = 1, ny
	READ(fgrid,*)dummy, yyv(j)
	yv(:,j) = yyv(j)
 END DO
 CLOSE(fgrid)

 !Reset Domain limits
 Xmin = MINVAL(xxv); Xmax = MAXVAL(xxv)
 Ymin = MINVAL(yyv); Ymax = MAXVAL(yyv)
 Lx = Xmax-Xmin; Ly = Ymax-Ymin

 !Create cell-centers
 DO j = 1, nyc
	DO i = 1, nxc
		xc(i,j) = 0.5D0*(xv(i,j) + xv(i+1,j))
		yc(i,j) = 0.5D0*(yv(i,j) + yv(i,j+1))
	END DO
 END DO

 !Create grid points with ghost points
 DO j = 1, nyg
	xvg(1,j)   = xxv(1)  - (xxv(2)  - xxv(1))
	xvg(nxg,j) = xxv(nx) + (xxv(nx) - xxv(nx-1))
 END DO
 DO i = 1, nxg
	yvg(i,1)   = yyv(1)  - (yyv(2)  - yyv(1))
	yvg(i,nyg) = yyv(ny) + (yyv(nx) - yyv(nx-1))
 END DO
 
 DO j = 1, ny
	yvg(:,j+1) = yyv(j)
 END DO
 DO i = 1, nx
	xvg(i+1,:) = xxv(i)
 END DO

 !Create cell-centers with ghost cells
 DO j = 1, nycg
	DO i = 1, nxcg
		xcg(i,j) = 0.5D0*(xvg(i,j) + xvg(i+1,j))
		ycg(i,j) = 0.5D0*(yvg(i,j) + yvg(i,j+1))
	END DO
 END DO

 WRITE(*,'(A)') ' DONE!'
 WRITE(*,'(A)') ''

 END SUBROUTINE GRID_NON_UNIFORM
 !----------------------------------------------------------------------------!


 SUBROUTINE CREATE_GRID_DATA_STRUC

 USE GLOBAL
 USE MESH2D
 USE IBM_PARAMS
 IMPLICIT NONE

 INTEGER :: i, j, ic, jc, iG, jG, icG, jcG, ICELL, ICELLG

 DOUBLE PRECISION :: dxew, dyns

 WRITE(*,'(A)',ADVANCE='NO') ' CREATING GRID DATA STRUCTURES ...'

 ibE = 0; ibW = 0; ibN = 0; ibS = 0; IC2NB = 0; IC2NBG = 0

 !Create associated data structures
 DO jc = 1, nyc
	DO ic = 1, nxc
		ICELL = (jc-1)*nxc + ic

		IC2CEL(ic,jc) = ICELL
		ICEL2C(ICELL,1) = ic; ICEL2C(ICELL,2) = jc;
 	END DO
 END DO
 DO jG = 1, nycg
	DO iG = 1, nxcg
		ICELLG = (jG-1)*nxcg + iG

		IC2CELG(iG,jG) = ICELLG
		ICELG2C(ICELLG,1) = iG; ICELG2C(ICELLG,2) = jG;
 	END DO
 END DO

 !Create neighbor data structures
 DO jc = 1, nyc
	DO ic = 1, nxc
		IF (ic .NE. nxc)	IC2NB(ic,jc,1) = IC2CEL(ic+1,jc)		!East   neighbor
		IF (ic .NE. 1)		IC2NB(ic,jc,2) = IC2CEL(ic-1,jc)		!West   neighbor
		IF (jc .NE. nyc)	IC2NB(ic,jc,3) = IC2CEL(ic,jc+1)		!North  neighbor
		IF (jc .NE. 1)		IC2NB(ic,jc,4) = IC2CEL(ic,jc-1)		!South  neighbor
	END DO
 END DO
 DO jc = 1, nyc
	jG = jc + 1
	DO ic = 1, nxc
		iG = ic + 1
		IC2NBG(iG,jG,1) = IC2CELG(iG+1,jG)		!East   neighbor
		IC2NBG(iG,jG,2) = IC2CELG(iG-1,jG)		!West   neighbor
		IC2NBG(iG,jG,3) = IC2CELG(iG,jG+1)		!North  neighbor
		IC2NBG(iG,jG,4) = IC2CELG(iG,jG-1)		!South  neighbor
	END DO
 END DO

 !Create inverse-cell-size data structures
 DO jc = 1, nyc
	DO ic = 1, nxc
		iG = ic+1; jG=jc+1
		!Inverse dimensions of cell
		dxcinv(ic,jc) = 1.D0/(xv(ic+1,jc) - xv(ic,jc))
		dycinv(ic,jc) = 1.D0/(yv(ic,jc+1) - yv(ic,jc))
		!Inverse step-size between cell-centers (non-compact stencil)
		dxcinv_NC(ic,jc) = 1.D0/(xcg(iG+1,jG) - xcg(iG-1,jG))
		dycinv_NC(ic,jc) = 1.D0/(ycg(iG,jG+1) - ycg(iG,jG-1))
	END DO
 END DO

 !Create inverse grid-spacing across face-centers
 DO jc = 1, nyc
	DO i = 1, nx
		icG = i; jcG = jc+1
		dxfinv(i,jc) = 1.D0/(xcg(icG+1,jcG)-xcg(icG,jcG))
 	END DO
 END DO
 DO j = 1, ny
	DO ic = 1, nxc
		icG = ic+1; jcG = j
		dyfinv(ic,j) = 1.D0/(ycg(icG,jcG+1)-ycg(icG,jcG))
	END DO
 END DO

 !Create inverse-square-distance data structures
 DO jc = 1, nyc
	DO ic = 1, nxc
		dxew = xv(ic+1,jc) - xv(ic,jc)
		dyns = yv(ic,jc+1) - yv(ic,jc)

		iG = ic + 1; jG = jc + 1
		dxe2inv(ic,jc) = 1.D0/dxew/( xcg(iG+1,jG) - xcg(iG,jG) )
		dxw2inv(ic,jc) = 1.D0/dxew/( xcg(iG,jG) - xcg(iG-1,jG) )
		dyn2inv(ic,jc) = 1.D0/dyns/( ycg(iG,jG+1) - ycg(iG,jG) )
		dys2inv(ic,jc) = 1.D0/dyns/( ycg(iG,jG) - ycg(iG,jG-1) )
	END DO
 END DO

 !Boundary Markers
 ibE(nxc,:) = 1; ibW(1,:) = 1
 ibN(:,nyc) = 1; ibS(:,1) = 1

 IF (iIntlBody) THEN
	DO jc = 1, nyc
		DO ic = 1, nxc
			icG = ic + 1; jcG = jc + 1
			IF (IBLANK(icG,jcG) .EQ. FLUID) THEN
				IF (IBLANK(icG+1,jcG) .EQ. SOLID) ibE(ic,jc) = 1
				IF (IBLANK(icG-1,jcG) .EQ. SOLID) ibW(ic,jc) = 1
				IF (IBLANK(icG,jcG+1) .EQ. SOLID) ibN(ic,jc) = 1
				IF (IBLANK(icG,jcG-1) .EQ. SOLID) ibS(ic,jc) = 1
			END IF
		END DO
 	END DO
 END IF

 WRITE(*,'(A)') ' DONE!'
 WRITE(*,'(A)') ''
 
 END SUBROUTINE CREATE_GRID_DATA_STRUC
 !----------------------------------------------------------------------------!


 SUBROUTINE DISP_GRID_STATS

 USE GLOBAL
 USE MESH2D
 IMPLICIT NONE

 WRITE(*,*) '------------------------------------------------------------------------'
 WRITE(*,*) 'GRID STATISTICS:'
 WRITE(*,*) ''
 WRITE(*,'(2(A,I3))') ' >> NUMBER OF INTERNAL GRID POINTS ALONG X & Y: ', nx, ',', ny
 WRITE(*,'(2(A,E11.4))') ' >> GRID SPACING: ', dx, ',', dy
 WRITE(*,*) '------------------------------------------------------------------------'
 WRITE(*,*) ''

 END SUBROUTINE DISP_GRID_STATS
 !----------------------------------------------------------------------------!
