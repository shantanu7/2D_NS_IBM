 SUBROUTINE FLUX_DIFFUSIVE

 USE GLOBAL
 USE MESH2D
 USE FLOW_PARAMS
 IMPLICIT NONE

 DOUBLE PRECISION :: dx2u, dy2u, dx2v, dy2v

 INTEGER :: ic, jc, icG, jcG

 DO jcG = 2, nycg-1
	DO icG = 2, nxcg-1

		ic = icG-1; jc = jcG-1

		!Diffusive Flux for x-momentum equation
		dx2u =	dxe2inv(ic,jc)*(uc(icG+1,jcG)-uc(icG,jcG)) - &
				dxw2inv(ic,jc)*(uc(icG,jcG)-uc(icG-1,jcG))
		dy2u =	dyn2inv(ic,jc)*(uc(icG,jcG+1)-uc(icG,jcG)) - &
				dys2inv(ic,jc)*(uc(icG,jcG)-uc(icG,jcG-1))
		FdifU(icG,jcG) = Reinv*(dx2u + dy2u)
		
		!Diffusive Flux for y-momentum equation
		dx2v =	dxe2inv(ic,jc)*(vc(icG+1,jcG)-vc(icG,jcG)) - &
				dxw2inv(ic,jc)*(vc(icG,jcG)-vc(icG-1,jcG))
		dy2v =	dyn2inv(ic,jc)*(vc(icG,jcG+1)-vc(icG,jcG)) - &
				dys2inv(ic,jc)*(vc(icG,jcG)-vc(icG,jcG-1))
		FdifV(icG,jcG) = Reinv*(dx2v + dy2v)

	END DO
 END DO

 END SUBROUTINE FLUX_DIFFUSIVE
 !----------------------------------------------------------------------------!
