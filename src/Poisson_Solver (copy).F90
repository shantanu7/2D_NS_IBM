 SUBROUTINE ALLOCATE_POISSON_MEMORY

 USE GLOBAL
 USE POISSON_PARAMS
 USE MESH2D
 IMPLICIT NONE

 ALLOCATE(Amat(NCELL,NCELL),phi(NCELL),rhs(NCELL))
 ALLOCATE(PBC(NCELL))

 END SUBROUTINE ALLOCATE_POISSON_MEMORY
 !----------------------------------------------------------------------------!


 SUBROUTINE INIT_POISSON

 USE GLOBAL
 USE POISSON_PARAMS
 USE MESH2D
 IMPLICIT NONE

 DOUBLE PRECISION :: xe, xw, yn, ys
 DOUBLE PRECISION :: dxew, dyns

 INTEGER :: i, j, ic, jc, ICELL, JCELL
 INTEGER icE, icW, icN, icS

 WRITE(*,'(A)',ADVANCE='NO') ' INITIALIZING POISSON DATA STRUCTURES ...'

 phi = 0.D0

 !Assemble coefficient and rhs for interior points
 DO jc = 1, nyc
	DO ic = 1, nxc
		ICELL = IC2CEL(ic,jc)
		icE = IC2NB(ic,jc,1)
		icW = IC2NB(ic,jc,2)
		icN = IC2NB(ic,jc,3)
		icS = IC2NB(ic,jc,4)

		Amat(ICELL,ICELL) = -1.D0*( dxe2inv(ic,jc)+dxw2inv(ic,jc) + &
								    dyn2inv(ic,jc)+dys2inv(ic,jc) )
		IF (ic .NE. nxc) Amat(ICELL,icE) 	=	dxe2inv(ic,jc)
		IF (ic .NE. 1)	 Amat(ICELL,icW)	=	dxw2inv(ic,jc)
		IF (jc .NE. nyc) Amat(ICELL,icN)	= 	dyn2inv(ic,jc)
		IF (jc .NE. 1)	 Amat(ICELL,icS)	=	dys2inv(ic,jc)

		!Commented here because RHS is now set in subroutine set_rhsp
!		rhs(ICELL) = f_Uvar(ic,jc)

	END DO
 END DO


 IF (MAXVAL(PBCType) .EQ. 0) THEN
	ic = INT(nxc/2)+1; jc = 1
	ICELL = IC2CEL(ic,jc)
	
	Amat(ICELL,:)		= 0.D0
	Amat(ICELL,ICELL)	= 1.D0
 END IF

 WRITE(*,'(A)') ' DONE!'
 PRINT*,''

 END SUBROUTINE INIT_POISSON
 !----------------------------------------------------------------------------!


 SUBROUTINE SET_PBC

 USE GLOBAL
 USE POISSON_PARAMS
 USE MESH2D
 IMPLICIT NONE

 DOUBLE PRECISION :: dxew, dyns
 DOUBLE PRECISION :: phiE, phiW, phiN, phiS
 DOUBLE PRECISION :: PBCD, PBCN

 INTEGER :: i, j, ic, jc, iG, jG, ICELL, JCELL
 INTEGER icE, icW, icN, icS

! WRITE(*,'(A)',ADVANCE='NO') ' SETTING PBC ...'

 PBC  = 0.D0
 phiE = PBCVal(1)
 phiW = PBCVal(2)
 phiN = PBCVal(3)
 phiS = PBCVal(4)

 !Assemble coefficient and rhs for interior points
 DO ICELL = 1, NCELL

	ic = ICEL2C(ICELL,1); jc = ICEL2C(ICELL,2)
	iG = ic + 1; jG = jc + 1
		
	icE = IC2NB(ic,jc,1)
	icW = IC2NB(ic,jc,2)
	icN = IC2NB(ic,jc,3)
	icS = IC2NB(ic,jc,4)

	!--------------------------------------------------------------!
	!Note:
	!In this implementation, the boundary condition vector (PBC) is 
	!reflected only in the initial computation. So the A matrix is 
	!modified to include only the Neumann BC. This component is sub-
	!tracted from the PBC vector.
	!--------------------------------------------------------------!

	!Calculate Dirichlet Boundary Condition
	PBCD =	-1.D0*( REAL(ibE(ic,jc)*(PBCType(1)),KIND=8)*dxe2inv(ic,jc)*phiE + & 
					REAL(ibW(ic,jc)*(PBCType(2)),KIND=8)*dxw2inv(ic,jc)*phiW + &
					REAL(ibN(ic,jc)*(PBCType(3)),KIND=8)*dyn2inv(ic,jc)*phiN + &
					REAL(ibS(ic,jc)*(PBCType(4)),KIND=8)*dys2inv(ic,jc)*phiS )

	PBC(ICELL) = PBCD

	!Add Neumann BC to A matrix
	Amat(ICELL,ICELL) = Amat(ICELL,ICELL) + &
				  ( REAL(ibE(ic,jc)*(1-PBCType(1)),KIND=8)*dxe2inv(ic,jc) + &
					REAL(ibW(ic,jc)*(1-PBCType(2)),KIND=8)*dxw2inv(ic,jc) + &
					REAL(ibN(ic,jc)*(1-PBCType(3)),KIND=8)*dyn2inv(ic,jc) + &
					REAL(ibS(ic,jc)*(1-PBCType(4)),KIND=8)*dys2inv(ic,jc) )
 END DO

 !Force one Dirichlet point
 IF (MAXVAL(PBCType) .EQ. 0) THEN
	ic = INT(nxc/2)+1; jc = 1
	iG = ic + 1; jG = jc + 1

	ICELL = IC2CEL(ic,jc)

	Amat(ICELL,:) = 0.D0
	Amat(ICELL,ICELL) = 1.D0

	PBC(ICELL) = 0.D0
	RHS(ICELL) = 1.D0
 END IF

! WRITE(*,'(A)') ' DONE!'
! PRINT*,''

 END SUBROUTINE SET_PBC
 !----------------------------------------------------------------------------!


 SUBROUTINE SET_PRESSURE_BC

 USE GLOBAL
 USE POISSON_PARAMS
 USE FLOW_PARAMS
 USE MESH2D
 IMPLICIT NONE

 DOUBLE PRECISION :: dxew, dyns
 DOUBLE PRECISION :: phiE, phiW, phiN, phiS
 DOUBLE PRECISION :: PBCD, PBCN

 INTEGER :: i, j, ic, jc, icG, jcG, ICELL, JCELL
 INTEGER icE, icW, icN, icS

 !Right Boundary
 DO jcG = 2, nycg-1
	pc(nxcg,jcG) = pc(nxcg-1,jcG)
 END DO 

 !Left Boundary -- Dirichlet
 DO jcG = 2, nycg-1
	pc(1,jcG) = 2.D0*PBCVal(2) - pc(2,jcG)
 END DO

 !Top Boundary
 DO icG = 2, nxcg-1
	pc(icG,nycg) = pc(icG,nycg-1)
 END DO

 !Bottom Boundary
 DO icG = 2, nxcg-1
	pc(icG,1) = pc(icG,2)
 END DO

 END SUBROUTINE SET_PRESSURE_BC
 !----------------------------------------------------------------------------!


 SUBROUTINE SOLVE_POISSON

 USE GLOBAL
 USE MESH2D
 USE POISSON_PARAMS
 USE FLOW_PARAMS
 IMPLICIT NONE

 INTEGER :: ic, jc, icG, jcG, ICELL, iter_hold

 iter_hold = PITERMAX

 CALL SET_PBC

 IF (TRIM(P_Solver) .EQ. 'JAC')  THEN
	CALL ITERATE_JACOBI(Amat, phi, (rhs+PBC), iter_hold, P_tol)
 ELSEIF (TRIM(P_Solver) .EQ. 'SOR') THEN
	CALL ITERATE_SOR(Amat, phi, (rhs+PBC), iter_hold, P_tol)
 ELSEIF (TRIM(P_Solver) .EQ. 'CG') THEN
	CALL ITERATE_CG(Amat, phi, (rhs+PBC), iter_hold, P_tol)
 END IF

 DO ICELL = 1, NCELL
	icG = ICEL2C(ICELL,1)+1; jcG = ICEL2C(ICELL,2)+1
	pc(icG,jcG) = PHI(ICELL)
 END DO

 END SUBROUTINE SOLVE_POISSON
 !----------------------------------------------------------------------------!


 SUBROUTINE SOLVE_PRESSURE

 USE GLOBAL
 USE MESH2D
 USE FLOW_PARAMS
 USE POISSON_PARAMS
 USE IBM_PARAMS
 IMPLICIT NONE

 INTEGER :: ic, jc, icG, jcG, iCELL, ICELLG
 INTEGER :: icE, icW, icN, icS
 INTEGER :: it

 DOUBLE PRECISION :: aP, aW, aE, aN, aS
 DOUBLE PRECISION :: source, tempX, tempY, DBC
 DOUBLE PRECISION :: ERROR, L2Norm, tolerance
 DOUBLE PRECISION :: CFL, CFLmax

 LOGICAL :: FLAG

 FLAG = .FALSE.

 it = 0; CFLmax = 0.D0

 DO WHILE ((.NOT. FLAG) .AND. it .LT. PITERMAX)

	ERROR = 0.D0
 
	!Solve PPE system
	DO jc = 1, nyc
		DO ic = 1, nxc
			icG = ic + 1; jcG = jc + 1
			ICELL = IC2CEL(ic,jc)

			ICELLG = IC2CELG(icG,jcG)
			icE = IC2NBG(icG,jcG,1)
			icW = IC2NBG(icG,jcG,2)
			icN = IC2NBG(icG,jcG,3)
			icS = IC2NBG(icG,jcG,4)

			aE = AmatPG(ICELLG,icE)
			aW = AmatPG(ICELLG,icW)
			aN = AmatPG(ICELLG,icN)
			aS = AmatPG(ICELLG,icS)
			aP = AmatPG(ICELLG,ICELLG)

			source = RHS(ICELL)

			pc(icG,jcG) = (1.D0-omega)*pc(icG,jcG) + &
					omega*(source - (aE*pc(icG+1,jcG) + aW*pc(icG-1,jcG) + &
									 aN*pc(icG,jcG+1) + aS*pc(icG,jcG-1)))/aP
	 	END DO
	END DO

	CALL SOLVE_PRESSURE_BC

	!Compute Error
	DO jc = 1, nyc
		DO ic = 1, nxc
			icG = ic + 1; jcG = jc + 1
			IF (IBLANK(icG,jcG) .EQ. FLUID) THEN
				ICELL = IC2CEL(ic,jc)

				aE = dxe2inv(ic,jc)
				aW = dxw2inv(ic,jc)
				aN = dyn2inv(ic,jc)
				aS = dys2inv(ic,jc)
				aP =  -1.D0*(aE+aW+aN+aS)

				ERROR = RHS(ICELL) - (aP*pc(icG,jcG) + &
				aE*pc(icG+1,jcG) + aW*pc(icG-1,jcG) + aN*pc(icG,jcG+1) + aS*pc(icG,jcG-1))
				L2NORM = L2NORM + ERROR**2
			END IF
		END DO
	END DO

	L2Norm = SQRT(L2Norm/REAL(nxc*nyc))
!	PRINT*,'ITERATION', it, 'ERROR', L2Norm
	IF (L2Norm .LT. P_Tol) THEN
		FLAG = .TRUE.
		WRITE(*,'(A,I4)')'     PRESSURE POISSON CONVERGED AFTER ITERATION: ', it
	END IF

	it = it + 1

 END DO
 WRITE(*,'(A,E)') '     ERROR: ', L2NORM

 END SUBROUTINE SOLVE_PRESSURE
 !----------------------------------------------------------------------------!


 SUBROUTINE DEALLOCATE_POISSON_MEMORY

 USE POISSON_PARAMS
 IMPLICIT NONE

 DEALLOCATE(Amat,phi,rhs,PBC) 

 END SUBROUTINE DEALLOCATE_POISSON_MEMORY

