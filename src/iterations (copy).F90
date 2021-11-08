 SUBROUTINE ITERATIONS

 USE GLOBAL
 USE MESH2D
 USE FLOW_PARAMS
 USE POISSON_PARAMS
 IMPLICIT NONE
 
 DOUBLE PRECISION :: tolerance = 1.0E-5, ERROR(2), ERR_SUM
 
 INTEGER :: ic, jc, icG, jcG, ICELL, ivar, k
 
 LOGICAL :: flag_u, flag_v

 pc = 0.
!! CALL FD_EULER

! !Pack velocity vector

! CALL ASSEMBLE_VEL_SYSTEM

! CALL SET_VBC_VECTOR

 DO ITER = 1, ITERMAX

!	ERROR = 1.E10; ERR_SUM = 1.E10

!	CALL FLUX_ADVECTIVE
!	CALL FLUX_DIFFUSIVE

!	DO ICELL = 1, NCELL
!		ic = ICEL2C(ICELL,1); jc = ICEL2C(ICELL,2)
!		icG = ic + 1; jcG = jc + 1
!		rhsV(ICELL,1) = uc(icG,jcG) + (-FconU(icG,jcG) + 0.5D0*FdifU(icG,jcG))*dt
!		rhsV(ICELL,2) = vc(icG,jcG) + (-FconV(icG,jcG) + 0.5D0*FdifV(icG,jcG))*dt

!		phiV(ICELL,1) = uc(icG,jcG); phiV(ICELL,2) = vc(icG,jcG)
!	END DO

!	DO WHILE (ERR_SUM .GT. tolerance)

!		DO ivar = 1, 2
!			CALL ITERATE_SOR(AmatV(:,:,ivar), phiV(:,ivar), & 
!							 rhsV(:,ivar)+VBC(:,ivar), 1, ERROR(ivar))
!		END DO
!		PRINT*,'ERROR', ERROR
!		ERR_SUM = 0.5*(ERROR(1) + ERROR(2))

!		DO ICELL = 1, NCELL
!			icG = ICEL2C(ICELL,1)+ 1; jcG = ICEL2C(ICELL,2) + 1
!			uc(icG,jcG) = phiV(ICELL,1); vc(icG,jcG) = phiV(ICELL,2)
!		END DO

!	END DO

!	

!	!Update intermediate velocity on face centers
!!	CALL INTERPOLATE_FC_VELOCITY

!!	!Set RHS for the Poisson Solver
!!	CALL SET_RHSP

!!	CALL SOLVE_POISSON

!!	CALL UPDATE_VELOCITY

!!	CALL SET_VBC_CC; CALL SET_VBC_FC

	CALL SOLVE_VELOCITY

	ctime = ctime + dt

	IF (MOD(ITER,TECFREQ) .EQ. 0) 	CALL WRITE_TEC_OUT_ASCII


 END DO

 END SUBROUTINE ITERATIONS
 !----------------------------------------------------------------------------!


 SUBROUTINE FD_EULER

 USE GLOBAL
 USE MESH2D
 USE FLOW_PARAMS
 USE POISSON_PARAMS
 IMPLICIT NONE

 INTEGER :: i, j, ic, jc, icG, jcG

 DO ITER = 1, ITERMAX

	PRINT*,''
	PRINT*,'>> ITER: ', ITER

	!Compute advective and diffusive fluxes along x- and y-
	CALL FLUX_ADVECTIVE
	CALL FLUX_DIFFUSIVE

	pc = 0.

	!Update intermediate velocity in the absence of pressure flux
	usr = uc + (-FconU + FdifU)*dt
	vsr = vc + (-FconV + FdifV)*dt
!	usr = uc -FconU*dt
!	vsr = vc -FconV*dt

	CALL SET_VBC_SR
	CALL SET_VBC_FC_SR

	!Update intermediate velocity on face centers
	CALL INTERPOLATE_FC_VELOCITY

	!Set RHS for the Poisson Solver
!	CALL SET_RHSP

!	CALL SOLVE_POISSON

	!Update velocity at cell-centers
	CALL UPDATE_VELOCITY

	CALL SET_VBC_CC; CALL SET_VBC_FC

	ctime = ctime + dt

	IF (MOD(ITER,TECFREQ) .EQ. 0) CALL WRITE_TEC_OUT_BINARY
 END DO

 PRINT*,''

 END SUBROUTINE FD_EULER
 !----------------------------------------------------------------------------!


 SUBROUTINE ASSEMBLE_VEL_SYSTEM

 USE GLOBAL
 USE MESH2D
 USE FLOW_PARAMS
 IMPLICIT NONE

 DOUBLE PRECISION :: xe, xw, yn, ys
 DOUBLE PRECISION :: dxew, dyns

 INTEGER :: i, j, ic, jc, ICELL, JCELL
 INTEGER icE, icW, icN, icS

 WRITE(*,'(A)',ADVANCE='NO') ' INITIALIZING VELOCITY DATA STRUCTURES ...'

 !Assemble coefficient and rhs for interior points
 DO jc = 1, nyc
	DO ic = 1, nxc
		ICELL = IC2CEL(ic,jc)
		icE = IC2NB(ic,jc,1)
		icW = IC2NB(ic,jc,2)
		icN = IC2NB(ic,jc,3)
		icS = IC2NB(ic,jc,4)


		IF (ic .NE. nxc) AmatV(ICELL,icE,:) =	-0.5D0*dt*Reinv*dxe2inv(ic,jc)
		IF (ic .NE. 1)	 AmatV(ICELL,icW,:) =	-0.5D0*dt*Reinv*dxw2inv(ic,jc)
		IF (jc .NE. nyc) AmatV(ICELL,icN,:) = 	-0.5D0*dt*Reinv*dyn2inv(ic,jc)
		IF (jc .NE. 1)	 AmatV(ICELL,icS,:) =	-0.5D0*dt*Reinv*dys2inv(ic,jc)

		AmatV(ICELL,ICELL,:) = 1.D0 + 0.5D0*dt*Reinv*( dxe2inv(ic,jc) + dxw2inv(ic,jc) + &
													   dyn2inv(ic,jc) + dys2inv(ic,jc) ) 
	END DO
 END DO

 WRITE(*,'(A)') ' DONE!'
 PRINT*,''

 END SUBROUTINE ASSEMBLE_VEL_SYSTEM
 !----------------------------------------------------------------------------!


 SUBROUTINE INTERPOLATE_FC_VELOCITY

 USE GLOBAL
 USE MESH2D
 USE FLOW_PARAMS
 USE POISSON_PARAMS
 IMPLICIT NONE
 
 INTEGER :: i, j, ic, jc, icG, jcG
 
 DO jc = 1, nyc
	DO i = 1, nx
		icG = i; jcG = jc + 1
		Ufsr(i,jc) = 0.5D0*(usr(icG+1,jcG)+usr(icG,jcG)) !!!!Need to interpolate non-uniform grid
 	END DO
 END DO
 DO j = 1, ny
	DO ic = 1, nxc
		icG = ic+1; jcG = j
		Vfsr(ic,j) = 0.5D0*(vsr(icG,jcG+1)+vsr(icG,jcG)) !!!!Need to interpolate non-uniform grid
 	END DO
 END DO

 END SUBROUTINE INTERPOLATE_FC_VELOCITY
 !----------------------------------------------------------------------------!


 SUBROUTINE SET_RHSP

 USE GLOBAL
 USE MESH2D
 USE FLOW_PARAMS
 USE POISSON_PARAMS
 IMPLICIT NONE

 INTEGER :: i, j, ic, jc, icG, jcG, ICELL

 DO jc = 1, nyc
	DO ic = 1, nxc
		icG = ic + 1; jcG = jc + 1
		ICELL = IC2CEL(ic,jc)

		RHS(ICELL) =  (	(Ufsr(ic+1,jc)-Ufsr(ic,jc))*dxcinv(ic,jc) + &
						(Vfsr(ic,jc+1)-Vfsr(ic,jc))*dycinv(ic,jc) )*dtinv

!		PHI(ICELL) = pc(icG,jcG)

!		PRINT*,'RHS', ic, jc, RHS(ICELL)
	END DO
 END DO

 END SUBROUTINE SET_RHSP
 !----------------------------------------------------------------------------!


 SUBROUTINE UPDATE_VELOCITY

 USE GLOBAL
 USE MESH2D
 USE FLOW_PARAMS
 USE POISSON_PARAMS
 IMPLICIT NONE

 INTEGER :: i, j, ic, jc, icG, jcG, ICELL

 !Update inner cell-center velocities
 !u^(n+1) = u* - dt*grad(p) -- non-compact
 DO jc = 1, nyc
	DO ic = 1, nxc
		icG = ic+1; jcG = jc+1
		uc(icG,jcG) = 	usr(icG,jcG) - &
						dt*dxcinv_NC(ic,jc)*(pc(icG+1,jcG)-pc(icG-1,jcG))
		vc(icG,jcG) = 	vsr(icG,jcG) - &
						dt*dycinv_NC(ic,jc)*(pc(icG,jcG+1)-pc(icG,jcG-1))
	END DO
 END DO

 !Update face-center velocities
 !Uf^(n+1) = Uf* - dt*grad(p) -- compact
 DO jc = 1, nyc
	DO i = 2, nx-1
		icG = i; jcG = jc+1
		Uf(i,jc) = Ufsr(i,jc) - dt*(pc(icG+1,jcG)-pc(icG,jcG))*dxfinv(i,jc)
 	END DO
 END DO
 DO j = 2, ny-1
	DO ic = 1, nxc
		icG = ic+1; jcG = j
		Vf(ic,j) = Vfsr(ic,j) - dt*(pc(icG,jcG+1)-pc(icG,jcG))*dyfinv(ic,j)
	END DO
 END DO

 END SUBROUTINE UPDATE_VELOCITY
 !----------------------------------------------------------------------------!


 SUBROUTINE SOLVE_VELOCITY

 USE GLOBAL
 USE MESH2D
 USE FLOW_PARAMS
 USE POISSON_PARAMS
 IMPLICIT NONE

 INTEGER :: ic, jc, icG, jcG, it

 DOUBLE PRECISION :: aP, aW, aE, aN, aS
 DOUBLE PRECISION :: source, tempX, tempY, DBC
 DOUBLE PRECISION :: ERROR, tolerance
 DOUBLE PRECISION :: CFL, CFLmax

 LOGICAL :: FLAG

 tolerance = 1.0E-8; FLAG = .FALSE.; pc = 0.

 !Set BC
 !Right Face:
 DO jcG = 2, nycg-1
	uc(nxcg,jcG) = uc(nxcg-1,jcG)
	vc(nxcg,jcG) = vc(nxcg-1,jcG)
 END DO

 !Left Face:
 DO jcG = 2, nycg-1
	DBC = 1.D0 - 4.D0*BCVal(2)*ycg(1,jcG)**2
	uc(1,jcG) = 2.D0*DBC - uc(2,jcG)
	vc(1,jcG) = -1.D0*vc(2,jcG)
 END DO

 !Top face:
 DO icG = 2, nxcg-1
	uc(icG,nycg) = 2.D0*BCVal(3) -uc(icG,nycg-1)
	vc(icG,nycg) = -1.D0*vc(icG,nycg-1)
 END DO

 !Bottom Face:
 DO icG = 2, nxcg-1
	uc(icG,1) = 2.D0*BCVal(4) -uc(icG,2)
	vc(icG,1) = -1.D0*vc(icG,2)
 END DO


 CALL FLUX_ADVECTIVE
 CALL FLUX_DIFFUSIVE

 it = 0; CFLmax = 0.D0

 DO WHILE ((.NOT. FLAG) .AND. it .LE. 100)

	ERROR = 0.D0
 
	!Solve CN system
	DO jc = 1, nyc
		DO ic = 1, nxc
			icG = ic + 1; jcG = jc + 1

			aE = -0.5D0*dt*Reinv*dxe2inv(ic,jc)
			aW = -0.5D0*dt*Reinv*dxw2inv(ic,jc)
			aN = -0.5D0*dt*Reinv*dyn2inv(ic,jc)
			aS = -0.5D0*dt*Reinv*dys2inv(ic,jc)
			aP = 1.D0 - (aE+aW+aN+aS) 

			tempX  = usr(icG,jcG); tempY = vsr(icG,jcG)

			source = uc(icG,jcG) + (-(1.5*FconU(icG,jcG)-0.5*FconU_old(icG,jcG)) + &
									  0.5*FdifU(icG,jcG))*dt

			usr(icG,jcG) = (source - (aE*usr(icG+1,jcG) + aW*usr(icG-1,jcG) + &
									  aN*usr(icG,jcG+1) + aS*usr(icG,jcG-1)))/aP

			source = vc(icG,jcG) + (-(1.5*FconV(icG,jcG)-0.5*FconV_old(icG,jcG)) + &
									  0.5*FdifV(icG,jcG))*dt

			vsr(icG,jcG) = (source - (aE*vsr(icG+1,jcG) + aW*vsr(icG-1,jcG) + &
									  aN*vsr(icG,jcG+1) + aS*vsr(icG,jcG-1)))/aP

			ERROR = ERROR + (usr(icG,jcG)-tempX)**2 + (vsr(icG,jcG)-tempY)**2
	 	END DO
	END DO

	CALL SOLVE_VELOCITY_BC_SR

	ERROR = SQRT(ERROR/REAL(nxc*nyc))
	PRINT*,'ITERATION', it, 'ERROR', ERROR
	IF (ERROR .LT. Tolerance) FLAG = .TRUE.

	it = it + 1
 END DO

 !Update intermediate velocity on face centers
 CALL INTERPOLATE_FC_VELOCITY

 !Set RHS for the Poisson Solver
 CALL SET_RHSP

 CALL SOLVE_PRESSURE

 CALL UPDATE_VELOCITY

 CALL SET_VBC_CC;  CALL SET_VBC_FC

 DO jcG = 2, nycg-1
	DO icG = 2, nxcg-1
		ic = icG-1; jc = jcG-1
		CFL = ABS(uc(icG,jcG))*dt*dxcinv(ic,jc)
		CFLmax = MAX(CFL,CFLmax)
	END DO
 END DO

 PRINT*,'CFLmax', CFLmax

 END SUBROUTINE SOLVE_VELOCITY
 !----------------------------------------------------------------------------!


 SUBROUTINE SOLVE_VELOCITY_BC_SR

 USE GLOBAL
 USE MESH2D
 USE FLOW_PARAMS
 USE POISSON_PARAMS
 IMPLICIT NONE

 INTEGER :: ic, jc, icG, jcG

 DOUBLE PRECISION :: DBC

 !Set BC
 !Right Face:
 DO jcG = 2, nycg-1
	usr(nxcg,jcG) = usr(nxcg-1,jcG)
	vsr(nxcg,jcG) = vsr(nxcg-1,jcG)
 END DO

 !Left Face:
 DO jcG = 2, nycg-1
	DBC = 1.D0 - 4.D0*BCVal(2)*ycg(1,jcG)**2
	usr(1,jcG) = 2.D0*DBC - usr(2,jcG)
	vsr(1,jcG) = - vsr(2,jcG)
 END DO

 !Top face:
 DO icG = 2, nxcg-1
	usr(icG,nycg) = 2.D0*BCVal(3) -usr(icG,nycg-1)
	vsr(icG,nycg) = -vsr(icG,nycg-1)
 END DO

 !Bottom Face:
 DO icG = 2, nxcg-1
	usr(icG,1) = 2.D0*BCVal(4) -usr(icG,2)
	vsr(icG,1) = -vsr(icG,2)
 END DO

 END SUBROUTINE SOLVE_VELOCITY_BC_SR
 !----------------------------------------------------------------------------!


 SUBROUTINE SOLVE_PRESSURE_BC

 USE GLOBAL
 USE MESH2D
 USE FLOW_PARAMS
 USE POISSON_PARAMS
 IMPLICIT NONE

 INTEGER :: ic, jc, icG, jcG

 DOUBLE PRECISION :: DBC

 !Set BC
 !Right Face:
 !Calculate pressure at outlet: Hagen-Poiseuille condition
 DO jcG = 2, nycg-1
 	DBC = PBCVal(1) + (xcg(nxcg,jcG)-Xmin)*Reinv*(1.D0-8.D0/Ly**2)
	pc(nxcg,jcG) = DBC
 END DO

 !Left Face:
 DO jcG = 2, nycg-1
	pc(1,jcG) = 2.D0*PBCVal(2) - pc(2,jcG)
 END DO

 !Top face:
 DO icG = 2, nxcg-1
	pc(icG,nycg) = pc(icG,nycg-1)
 END DO

 !Bottom Face:
 DO icG = 2, nxcg-1
	pc(icG,1) = pc(icG,2)
 END DO

 END SUBROUTINE SOLVE_PRESSURE_BC
 !----------------------------------------------------------------------------!
