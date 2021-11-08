 SUBROUTINE FLUX_ADVECTIVE

 USE GLOBAL
 USE MESH2D
 USE FLOW_PARAMS
 IMPLICIT NONE

 DOUBLE PRECISION :: ue, uw, un, us
 DOUBLE PRECISION :: ve, vw, vn, vs
 DOUBLE PRECISION :: duUfx, duVfy, dvVfy, dvUfx


 INTEGER :: ic, jc, icG, jcG

 FconU_old = FconU; FconV_old = FconV

 DO jcG = 2, nycg-1
	DO icG = 2, nxcg-1
		ic = icG-1; jc = jcG-1

		!!!!Need to interpolate non-uniform grid
		ue	= 0.5*(uc(icG,jcG) + uc(icG+1,jcG))
		uw	= 0.5*(uc(icG,jcG) + uc(icG-1,jcG))
		un	= 0.5*(uc(icG,jcG) + uc(icG,jcG+1))
		us	= 0.5*(uc(icG,jcG) + uc(icG,jcG-1))

		ve	= 0.5*(vc(icG,jcG) + vc(icG+1,jcG))
		vw	= 0.5*(vc(icG,jcG) + vc(icG-1,jcG))
		vn	= 0.5*(vc(icG,jcG) + vc(icG,jcG+1))
		vs	= 0.5*(vc(icG,jcG) + vc(icG,jcG-1))

		duUfx = (ue*Uf(ic+1,jc)-uw*Uf(ic,jc))*dxcinv(ic,jc)
		duVfy = (un*Vf(ic,jc+1)-us*Vf(ic,jc))*dycinv(ic,jc)
		dvUfx = (ve*Uf(ic+1,jc)-vw*Uf(ic,jc))*dxcinv(ic,jc)
		dvVfy = (vn*Vf(ic,jc+1)-vs*Vf(ic,jc))*dycinv(ic,jc)

		!Convective Flux for x-momentum equation
		FconU(icG,jcG) = duUfx + duVfy

		!Convective Flux for y-momentum equation
		FconV(icG,jcG) = dvUfx + dvVfy
	END DO
 END DO

 IF (ITER .EQ. 1) THEN
	FconU_old = FconU; FconV_old = FconV
 END IF

 END SUBROUTINE FLUX_ADVECTIVE
 !----------------------------------------------------------------------------!
