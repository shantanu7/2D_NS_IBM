 SUBROUTINE ITERATIONS

 USE GLOBAL
 USE MESH2D
 USE FLOW_PARAMS
 USE POISSON_PARAMS
 IMPLICIT NONE
 
 DOUBLE PRECISION :: CFL, CFLmax
 
 INTEGER :: ic, jc, icG, jcG, ICELL, ivar, k
 
 LOGICAL :: flag_u, flag_v

 CALL ASSEMBLE_AMATG

 CFLmax = 0.D0

 WRITE(*,'(A)') ' STARTING ITERATIONS ...'

 DO ITER = 1, ITERMAX

	PRINT*,'!------------------------------------------------------------------!'
	WRITE(*,'(A,I5)')' >> ITER: ', ITER
	PRINT*,''

	CALL SOLVE_ADVECTION_DIFFUSION

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

	WRITE(*,'(A,E)')'     CFLmax:', CFLmax

	IF (CFLmax .GT. 1.D0) THEN
		WRITE(*,'(A,F10.7,A)') ' ERROR: CFLmax (', CFLmax ,') GREATER THAN 1.0!!!'
		CALL WRITE_TEC_OUT_ASCII
		WRITE(*,'(A)') ' ABORTING!!!'
		STOP
	END IF

	ctime = ctime + dt

	CALL WRITE_PROBE_DATA

	PRINT*,'!------------------------------------------------------------------!'

	IF (MOD(ITER,TECFREQ) .EQ. 0) 	CALL WRITE_TEC_OUT_ASCII

	PRINT*,''
	PRINT*,''
	PRINT*,''

 END DO

 CALL CLOSE_PROBE_FILES

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
	WRITE(*,'(A,I5)')' >> ITER: ', ITER

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

 WRITE(*,'(A)',ADVANCE='NO') ' INITIALIZING GLOBAL MATRICES...'

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


 SUBROUTINE ASSEMBLE_AMATG

 USE GLOBAL
 USE MESH2D
 USE FLOW_PARAMS
 USE IBM_PARAMS
 IMPLICIT NONE

 DOUBLE PRECISION :: xe, xw, yn, ys
 DOUBLE PRECISION :: dxew, dyns

 INTEGER :: i, j, ic, jc, icG, jcG, ICELLG
 INTEGER :: icE, icW, icN, icS
 INTEGER :: iSolid, iGC, iNB

 !Assemble part of Amat for fluid cells
 DO jc = 1, nyc
	DO ic = 1, nxc
		icG = ic + 1; jcG = jc + 1

		IF (IBLANK(icG,jcG) .EQ. FLUID) THEN
			ICELLG = IC2CELG(icG,jcG)
			icE = IC2NBG(icG,jcG,1)
			icW = IC2NBG(icG,jcG,2)
			icN = IC2NBG(icG,jcG,3)
			icS = IC2NBG(icG,jcG,4)

			!Assemble velocity matrix
			AmatVG(ICELLG,icE) = -0.5D0*dt*Reinv*dxe2inv(ic,jc)
			AmatVG(ICELLG,icW) = -0.5D0*dt*Reinv*dxw2inv(ic,jc)
			AmatVG(ICELLG,icN) = -0.5D0*dt*Reinv*dyn2inv(ic,jc)
			AmatVG(ICELLG,icS) = -0.5D0*dt*Reinv*dys2inv(ic,jc)

			AmatVG(ICELLG,ICELLG) = 1.D0 + 0.5D0*dt*Reinv* &
						( dxe2inv(ic,jc) + dxw2inv(ic,jc) + &
						  dyn2inv(ic,jc) + dys2inv(ic,jc) )

			!Assemble Pressure matrix
			AmatPG(ICELLG,icE) = dxe2inv(ic,jc)
			AmatPG(ICELLG,icW) = dxw2inv(ic,jc)
			AmatPG(ICELLG,icN) = dyn2inv(ic,jc)
			AmatPG(ICELLG,icS) = dys2inv(ic,jc)

			AmatPG(ICELLG,ICELLG) = -1.D0*( dxe2inv(ic,jc) + dxw2inv(ic,jc) + &
											dyn2inv(ic,jc) + dys2inv(ic,jc) )
		END IF

	END DO
 END DO

 !Now treat the Solid cells
 DO iSolid = 1, NSolid
	icG = IS2C(iSolid,1); jcG = IS2C(iSolid,2)
	ICELLG = IC2CELG(icG,jcG)

	!Clear previously assigned values
	AmatVG(ICELLG,:) = 0.D0
	AmatVG(ICELLG,ICELLG) = 1.D0	
	AmatPG(ICELLG,:) = 0.D0
	AmatPG(ICELLG,ICELLG) = 1.D0
 END DO

 !Add boundary condition for ghost cells
 DO iGC = 1, NGC
	icG = IGC2C(IGC,1); jcG = IGC2C(IGC,2)

	ICELLG	= IC2CELG(icG,jcG)
	iNB		= IC2CELG(IGC2FNB(iGC,1),IGC2FNB(iGC,2))

	!Add no slip BC
	AmatVG(ICELLG,iNB) = 1.D0
	!Add homogenous Neumann BC
	AmatPG(ICELLG,iNB) = -1.D0
 END DO

 END SUBROUTINE ASSEMBLE_AMATG 
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


 SUBROUTINE SET_RHSV

 USE GLOBAL
 USE MESH2D
 USE FLOW_PARAMS
 USE POISSON_PARAMS
 USE IBM_PARAMS
 IMPLICIT NONE

 INTEGER :: i, j, ic, jc, icG, jcG, ICELL

 RHSV = 0.D0

 DO jc = 1, nyc
	DO ic = 1, nxc
		icG = ic + 1; jcG = jc + 1

		IF (IBLANK(icG,jcG) .EQ. FLUID) THEN
			ICELL = IC2CEL(ic,jc)
			RHSV(ICELL,1) = uc(icG,jcG) + &
			( -(1.5*FconU(icG,jcG) - 0.5*FconU_old(icG,jcG)) + &
			    0.5*FdifU(icG,jcG) )*dt

			RHSV(ICELL,2) = vc(icG,jcG) + &
			( -(1.5*FconV(icG,jcG) - 0.5*FconV_old(icG,jcG)) + &
			    0.5*FdifV(icG,jcG) )*dt
		END IF
	END DO
 END DO

 END SUBROUTINE SET_RHSV
 !----------------------------------------------------------------------------!


 SUBROUTINE SET_RHSP

 USE GLOBAL
 USE MESH2D
 USE FLOW_PARAMS
 USE POISSON_PARAMS
 USE IBM_PARAMS
 IMPLICIT NONE

 INTEGER :: i, j, ic, jc, icG, jcG, ICELL

 RHS = 0.D0

 DO jc = 1, nyc
	DO ic = 1, nxc
		icG = ic + 1; jcG = jc + 1

		IF (IBLANK(icG,jcG) .EQ. FLUID) THEN
			ICELL = IC2CEL(ic,jc)

			RHS(ICELL) =  (	(Ufsr(ic+1,jc)-Ufsr(ic,jc))*dxcinv(ic,jc) + &
							(Vfsr(ic,jc+1)-Vfsr(ic,jc))*dycinv(ic,jc) )*dtinv
		END IF
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


 SUBROUTINE SOLVE_ADVECTION_DIFFUSION

 USE GLOBAL
 USE MESH2D
 USE FLOW_PARAMS
 USE POISSON_PARAMS
 IMPLICIT NONE

 INTEGER :: ic, jc, icG, jcG, ICELL, ICELLG
 INTEGER :: icE, icW, icN, icS
 INTEGER :: it

 DOUBLE PRECISION :: aP, aW, aE, aN, aS
 DOUBLE PRECISION :: source, tempX, tempY, DBC
 DOUBLE PRECISION :: ERROR, tolerance

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
	DBC = BCVal(2)*(1.D0 - 4.D0*(ycg(1,jcG)/Ly)**2)
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

 CALL SET_RHSV

 it = 0

 DO WHILE ((.NOT. FLAG) .AND. it .LE. 100)

	ERROR = 0.D0
 
	!Solve CN system
	DO jc = 1, nyc
		DO ic = 1, nxc
			icG = ic + 1; jcG = jc + 1
			
			ICELL  = IC2CEL(ic,jc)
			ICELLG = IC2CELG(icG,jcG)
			icE = IC2NBG(icG,jcG,1)
			icW = IC2NBG(icG,jcG,2)
			icN = IC2NBG(icG,jcG,3)
			icS = IC2NBG(icG,jcG,4)

			aE = AmatVG(ICELLG,icE)
			aW = AmatVG(ICELLG,icW)
			aN = AmatVG(ICELLG,icN)
			aS = AmatVG(ICELLG,icS)
			aP = AmatVG(ICELLG,ICELLG)

			tempX  = usr(icG,jcG); tempY = vsr(icG,jcG)

			!Solve x-Advection-Diffusion equation
			usr(icG,jcG) = (RHSV(ICELL,1) - (aE*usr(icG+1,jcG) + aW*usr(icG-1,jcG) + &
									  		 aN*usr(icG,jcG+1) + aS*usr(icG,jcG-1)))/aP

			!Solve y-Advection-Diffusion equation
			vsr(icG,jcG) = (RHSV(ICELL,2) - (aE*vsr(icG+1,jcG) + aW*vsr(icG-1,jcG) + &
									  		 aN*vsr(icG,jcG+1) + aS*vsr(icG,jcG-1)))/aP

			!Check combined error
			ERROR = ERROR + (usr(icG,jcG)-tempX)**2 + (vsr(icG,jcG)-tempY)**2
	 	END DO
	END DO

	CALL SOLVE_VELOCITY_BC_SR

	ERROR = SQRT(ERROR/REAL(nxc*nyc))

	IF (ERROR .LT. Tolerance) THEN
		FLAG = .TRUE.
		WRITE(*,'(A,I4)')'     ADVECTION-DIFFUSION CONVERGED AFTER ITERATION: ', it
		WRITE(*,'(A,E)') '     ERROR: ', ERROR
	END IF

	it = it + 1
 END DO

 END SUBROUTINE SOLVE_ADVECTION_DIFFUSION
 !----------------------------------------------------------------------------!


 SUBROUTINE SOLVE_VELOCITY_BC_SR

 USE GLOBAL
 USE MESH2D
 USE FLOW_PARAMS
 USE POISSON_PARAMS
 USE IBM_PARAMS, ONLY: D_cyl
 IMPLICIT NONE

 INTEGER :: ic, jc, icG, jcG

 DOUBLE PRECISION :: DBC, Umax

 !Set BC
 !Right Face:
 DO jcG = 2, nycg-1
	usr(nxcg,jcG) = usr(nxcg-1,jcG)
	vsr(nxcg,jcG) = vsr(nxcg-1,jcG)
 END DO

 !Left Face:
 DO jcG = 2, nycg-1
	DBC = BCVal(2)*(1.D0 - 4.D0*(ycg(1,jcG)/Ly)**2)
	usr(1,jcG) = 2.D0*DBC - usr(2,jcG)
	vsr(1,jcG) = -1.D0*vsr(2,jcG)
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
	pc(nxcg,jcG) = 	2.*PBCVal(1) - pc(nxcg-1,jcG)
 END DO

 !Left Face:
 DO jcG = 2, nycg-1
!	pc(1,jcG) = 2.*PBCVal(2) - pc(2,jcG)
	pc(1,jcG) = pc(2,jcG)
 END DO

! pc(1,INT(nycg/2)) = 2.*PBCVal(2) - pc(2,INT(nycg/2))

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
