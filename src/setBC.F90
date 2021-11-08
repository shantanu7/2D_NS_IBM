 SUBROUTINE SET_VBC_CC

 USE GLOBAL
 USE MESH2D
 USE FLOW_PARAMS

 DOUBLE PRECISION :: dh, a, b, DBC, NBC, test

 INTEGER :: i, j, ic, jc, icG, jcG

 !--------------
 !Right Boundary
 !--------------
 IF (BCType(1) .EQ. INLET) THEN
	DO jcG = 2, nycg-1
		j	= jcG-1
		dh	= xcg(nxcg,jcG)-xcg(nxcg-1,jcG)
		a	= xv(nx,j)-xcg(nxcg-1,jcG); b = xcg(nxcg,jcG) - xv(nx,j)
		DBC = (BCVal(1)*dh - b*uc(nxcg-1,jcG))/a
		uc(nxcg,jcG) = DBC
		vc(nxcg,jcG) = 0.D0
	END DO
 ELSEIF (BCType(1) .EQ. OUTLET) THEN
	DO jcG = 2, nycg-1
		j	= jcG-1
!		dh	= xcg(nxcg,jcG)-xcg(nxcg-1,jcG)
!		a	= xv(nx,j)-xcg(nxcg-1,jcG); b = xcg(nxcg,jcG) - xv(nx,j)
!		NBC = (BCVal(1)*dh + uc(nxcg-1,jcG))
		uc(nxcg,jcG) = uc(nxcg-1,jcG)
		vc(nxcg,jcG) = vc(nxcg-1,jcG)
	END DO
 ELSEIF (BCType(1) .EQ. NOSLIP) THEN
	DO jcG = 2, nycg-1
		j	= jcG-1
		dh	= xcg(nxcg,jcG)-xcg(nxcg-1,jcG)
		a	= xv(nx,j)-xcg(nxcg-1,jcG); b = xcg(nxcg,jcG) - xv(nx,j)
		DBC = (BCVal(1)*dh - b*vc(nxcg-1,jcG))/a
		uc(nxcg,jcG) = -1.D0*uc(nxcg-1,jcG)		!Zero-penetration on normal velocity
		vc(nxcg,jcG) = DBC						!No slip applies tangential velocity
	END DO
 END IF

 !-------------
 !Left Boundary 
 !-------------
 IF (BCType(2) .EQ. INLET) THEN
	DO jcG = 2, nycg-1
		j	= jcG-1
		dh	= xcg(2,jcG)-xcg(1,jcG)
		a	= xcg(2,jcG)-xv(1,j); b = xv(1,j)-xcg(1,jcG)

		test= BCVal(2)*(1.D0 - 4.D0*(ycg(1,jcG)/Ly)**2)
		uc(1,jcG) = 2.*test - uc(2,jcG)
		vc(1,jcG) = -vc(2,jcG)
	END DO
 ELSEIF (BCType(2) .EQ. OUTLET) THEN
	DO jcG = 2, nycg-1
		j	= jcG-1
		dh	= xcg(2,jcG)-xcg(1,jcG)
		a	= xcg(2,jcG)-xv(1,j); b = xv(1,j)-xcg(1,jcG)
		NBC = (BCVal(2)*dh + uc(2,jcG))
		uc(1,jcG) = NBC
		vc(1,jcG) = 0.D0
	END DO
 ELSEIF (BCType(2) .EQ. NOSLIP) THEN
	DO jcG = 2, nycg-1
		j	= jcG-1
		dh	= xcg(2,jcG)-xcg(1,jcG)
		a	= xcg(2,jcG)-xv(1,j); b = xv(1,j)-xcg(1,jcG)
		DBC = (BCVal(2)*dh - b*vc(2,jcG))/a
		uc(1,jcG) = -1.D0*uc(2,jcG)
		vc(1,jcG) = DBC
	END DO
 END IF

 !------------
 !Top Boundary
 !------------
 IF (BCType(3) .EQ. INLET) THEN
	DO icG = 2, nxcg-1
		i	= icG-1
		dh	= ycg(icG,nycg)-ycg(icG,nycg-1)
		a	= yv(i,ny)-ycg(icG,nycg-1); b = ycg(icG,nycg) - yv(i,ny)
		DBC = (BCVal(3)*dh - b*vc(icG,nyc-1))/a
		uc(icG,nycg) = 0.D0
		vc(icG,nycg) = DBC
	END DO
 ELSEIF (BCType(3) .EQ. OUTLET) THEN
	DO icG = 2, nxcg-1
		i	= icG-1
		dh	= ycg(icG,nycg)-ycg(icG,nycg-1)
		a	= yv(i,ny)-ycg(icG,nycg-1); b = ycg(icG,nycg) - yv(i,ny)
		NBC = (BCVal(3)*dh + vc(icG,nyc-1))
		uc(icG,nycg) = 0.D0
		vc(icG,nycg) = NBC
	END DO
 ELSEIF (BCType(3) .EQ. NOSLIP) THEN
	DO icG = 2, nxcg-1
		i	= icG-1
!		dh	= ycg(icG,nycg)-ycg(icG,nycg-1)
!		a	= yv(i,ny)-ycg(icG,nycg-1); b = ycg(icG,nycg) - yv(i,ny)
!		DBC = (BCVal(3)*dh - b*uc(icG,nyc-1))/a
		uc(icG,nycg) = 2.*BCVal(3) - 1.D0*uc(icG,nycg-1)
		vc(icG,nycg) = -1.D0*vc(icG,nycg-1)
	END DO
 END IF

 !---------------
 !Bottom Boundary
 !---------------
 IF (BCType(4) .EQ. INLET) THEN
	DO icG = 2, nxcg-1
		i	= icG-1
		dh	= ycg(icG,2)-ycg(icG,1)
		a	= ycg(icG,2)-yv(i,1); b = yv(i,1)-ycg(icG,1) 
		DBC = (BCVal(4)*dh - b*vc(icG,2))/a
		uc(icG,1) = 0.D0
		vc(icG,1) = DBC
	END DO
 ELSEIF (BCType(4) .EQ. OUTLET) THEN
	DO icG = 2, nxcg-1
		i	= icG-1
		dh	= ycg(icG,2)-ycg(icG,1)
		a	= ycg(icG,2)-yv(i,1); b = yv(i,1)-ycg(icG,1) 
		NBC = (BCVal(4)*dh + vc(icG,2))
		uc(icG,1) = 0.D0
		vc(icG,1) = NBC
	END DO
 ELSEIF (BCType(4) .EQ. NOSLIP) THEN
	DO icG = 2, nxcg-1
		i	= icG-1
!		dh	= ycg(icG,2)-ycg(icG,1)
!		a	= ycg(icG,2)-yv(i,1); b = yv(i,1)-ycg(icG,1) 
!		DBC = (BCVal(4)*dh - b*uc(icG,2))/a
		uc(icG,1) = 2.*BCVal(4) - 1.D0*uc(icG,2)
		vc(icG,1) = -1.D0*vc(icG,2)
	END DO
 END IF

 END SUBROUTINE SET_VBC_CC
 !----------------------------------------------------------------------------!


 SUBROUTINE SET_VBC_SR

 USE GLOBAL
 USE MESH2D
 USE FLOW_PARAMS

 DOUBLE PRECISION :: dh, a, b, DBC, NBC, test

 INTEGER :: i, j, ic, jc, icG, jcG

 !--------------
 !Right Boundary
 !--------------
 IF (BCType(1) .EQ. INLET) THEN
	DO jcG = 2, nycg-1
		j	= jcG-1
		dh	= xcg(nxcg,jcG)-xcg(nxcg-1,jcG)
		a	= xv(nx,j)-xcg(nxcg-1,jcG); b = xcg(nxcg,jcG) - xv(nx,j)
		DBC = (BCVal(1)*dh - b*usr(nxcg-1,jcG))/a
		usr(nxcg,jcG) = DBC
		vsr(nxcg,jcG) = 0.D0
	END DO
 ELSEIF (BCType(1) .EQ. OUTLET) THEN
	DO jcG = 2, nycg-1
!		j	= jcG-1
!		dh	= xcg(nxcg,jcG)-xcg(nxcg-1,jcG)
!		a	= xv(nx,j)-xcg(nxcg-1,jcG); b = xcg(nxcg,jcG) - xv(nx,j)
!		NBC = (BCVal(1)*dh + usr(nxcg-1,jcG))
!		usr(nxcg,jcG) = NBC
!		vsr(nxcg,jcG) = 0.D0
		usr(nxcg,jcG) = usr(nxcg-1,jcG)
		vsr(nxcg,jcG) = vsr(nxcg-1,jcG)
	END DO
 ELSEIF (BCType(1) .EQ. NOSLIP) THEN
	DO jcG = 2, nycg-1
		j	= jcG-1
		dh	= xcg(nxcg,jcG)-xcg(nxcg-1,jcG)
		a	= xv(nx,j)-xcg(nxcg-1,jcG); b = xcg(nxcg,jcG) - xv(nx,j)
		DBC = (BCVal(1)*dh - b*vsr(nxcg-1,jcG))/a
		usr(nxcg,jcG) = -1.D0*usr(nxcg-1,jcG)		!Zero-penetration on normal velocity
		vsr(nxcg,jcG) = DBC						!No slip applies tangential velocity
	END DO
 END IF

 !-------------
 !Left Boundary 
 !-------------
 IF (BCType(2) .EQ. INLET) THEN
	DO jcG = 2, nycg-1
!		j	= jcG-1
!		dh	= xcg(2,jcG)-xcg(1,jcG)
!		a	= xcg(2,jcG)-xv(1,j); b = xv(1,j)-xcg(1,jcG)
!		test= 1. - 4.*BCVal(2)*ycg(1,jcG)**2
!		DBC = (test*dh - b*usr(2,jcG))/a
!!		DBC = (BCVal(2)*dh - b*usr(2,jcG))/a
!		usr(1,jcG) = DBC
!		vsr(1,jcG) = 0.D0
		DBC = BCVal(2)*(1.D0 - 4.D0*(ycg(1,jcG)/Ly)**2)
		usr(1,jcG) = 2.D0*DBC - usr(2,jcG)
		vsr(1,jcG) = - 1.D0*vsr(2,jcG)
	END DO
 ELSEIF (BCType(2) .EQ. OUTLET) THEN
	DO jcG = 2, nycg-1
		j	= jcG-1
		dh	= xcg(2,jcG)-xcg(1,jcG)
		a	= xcg(2,jcG)-xv(1,j); b = xv(1,j)-xcg(1,jcG)
		NBC = (BCVal(2)*dh + usr(2,jcG))
		usr(1,jcG) = NBC
		vsr(1,jcG) = 0.D0
	END DO
 ELSEIF (BCType(2) .EQ. NOSLIP) THEN
	DO jcG = 2, nycg-1
		j	= jcG-1
		dh	= xcg(2,jcG)-xcg(1,jcG)
		a	= xcg(2,jcG)-xv(1,j); b = xv(1,j)-xcg(1,jcG)
		DBC = (BCVal(2)*dh - b*vsr(2,jcG))/a
		usr(1,jcG) = -1.D0*usr(2,jcG)
		vsr(1,jcG) = DBC
	END DO
 END IF

 !------------
 !Top Boundary
 !------------
 IF (BCType(3) .EQ. INLET) THEN
	DO icG = 2, nxcg-1
		i	= icG-1
		dh	= ycg(icG,nycg)-ycg(icG,nycg-1)
		a	= yv(i,ny)-ycg(icG,nycg-1); b = ycg(icG,nycg) - yv(i,ny)
		DBC = (BCVal(3)*dh - b*vsr(icG,nyc-1))/a
		usr(icG,nycg) = 0.D0
		vsr(icG,nycg) = DBC
	END DO
 ELSEIF (BCType(3) .EQ. OUTLET) THEN
	DO icG = 2, nxcg-1
		i	= icG-1
		dh	= ycg(icG,nycg)-ycg(icG,nycg-1)
		a	= yv(i,ny)-ycg(icG,nycg-1); b = ycg(icG,nycg) - yv(i,ny)
		NBC = (BCVal(3)*dh + vsr(icG,nyc-1))
		usr(icG,nycg) = 0.D0
		vsr(icG,nycg) = NBC
	END DO
 ELSEIF (BCType(3) .EQ. NOSLIP) THEN
	DO icG = 2, nxcg-1
!		i	= icG-1
!		dh	= ycg(icG,nycg)-ycg(icG,nycg-1)
!		a	= yv(i,ny)-ycg(icG,nycg-1); b = ycg(icG,nycg) - yv(i,ny)
!		DBC = (BCVal(3)*dh - b*usr(icG,nyc-1))/a
!		usr(icG,nycg) = DBC
!		vsr(icG,nycg) = -1.D0*vsr(icG,nycg-1)
		usr(icG,nycg) = 2.D0*BCVal(3) - 1.D0*usr(icG,nycg-1)
		vsr(icG,nycg) = -1.D0*vsr(icG,nycg-1)
	END DO
 END IF

 !---------------
 !Bottom Boundary
 !---------------
 IF (BCType(4) .EQ. INLET) THEN
	DO icG = 2, nxcg-1
		i	= icG-1
		dh	= ycg(icG,2)-ycg(icG,1)
		a	= ycg(icG,2)-yv(i,1); b = yv(i,1)-ycg(icG,1) 
		DBC = (BCVal(4)*dh - b*vsr(icG,2))/a
		usr(icG,1) = 0.D0
		vsr(icG,1) = DBC
	END DO
 ELSEIF (BCType(4) .EQ. OUTLET) THEN
	DO icG = 2, nxcg-1
		i	= icG-1
		dh	= ycg(icG,2)-ycg(icG,1)
		a	= ycg(icG,2)-yv(i,1); b = yv(i,1)-ycg(icG,1) 
		NBC = (BCVal(4)*dh + vsr(icG,2))
		usr(icG,1) = 0.D0
		vsr(icG,1) = NBC
	END DO
 ELSEIF (BCType(4) .EQ. NOSLIP) THEN
	DO icG = 2, nxcg-1
!		i	= icG-1
!		dh	= ycg(icG,2)-ycg(icG,1)
!		a	= ycg(icG,2)-yv(i,1); b = yv(i,1)-ycg(icG,1) 
!		DBC = (BCVal(4)*dh - b*usr(icG,2))/a
!		usr(icG,1) = DBC
!		vsr(icG,1) = -1.D0*vsr(icG,2)
		usr(icG,1) = 2.D0*BCVal(4) -1.D0*usr(icG,2)
		vsr(icG,1) = -1.D0*vsr(icG,2)
	END DO
 END IF

 END SUBROUTINE SET_VBC_SR
 !----------------------------------------------------------------------------!


 SUBROUTINE SET_VBC_FC

 USE GLOBAL
 USE MESH2D
 USE FLOW_PARAMS

 DOUBLE PRECISION :: dh, a, b, DBC, NBC

 INTEGER :: i, j, ic, jc, icG, jcG

 !--------------
 !Right Boundary
 !--------------
 IF (BCType(1) .EQ. INLET) THEN
	Uf(nx,:) = BCVal(1)
 ELSEIF (BCType(1) .EQ. OUTLET) THEN
	DO jc = 1, nyc
		jcG = jc+1; j = jc
		dh  = xv(nx,j) - xcg(nxcg-1,jcG)
		Uf(nx,jc) = BCVal(1)*dh + uc(nxcg-1,jcG)
	END DO
 ELSEIF (BCType(1) .EQ. NOSLIP) THEN
	Uf(nx,:) = 0.D0
 END IF

 !-------------
 !Left Boundary
 !-------------
 IF (BCType(2) .EQ. INLET) THEN
	Uf(1,:) = BCVal(2)
	DO jc = 1, nyc
		Uf(1,jc) = BCVal(2)*(1.D0 - 4.D0*(ycg(1,jcG)/Ly)**2)
 	END DO
 ELSEIF (BCType(2) .EQ. OUTLET) THEN
	DO jc = 1, nyc
		jcG = jc+1; j = jc
		dh  = xcg(2,jcG) - xv(1,j)
		Uf(1,jc) = BCVal(2)*dh + uc(2,jcG)
	END DO
 ELSEIF (BCType(2) .EQ. NOSLIP) THEN
	Uf(1,:) = 0.D0
 END IF

 !------------
 !Top Boundary
 !------------
 IF (BCType(3) .EQ. INLET) THEN
	Vf(:,ny) = BCVal(3)
 ELSEIF (BCType(3) .EQ. OUTLET) THEN
	DO ic = 1, nxc
		icG = ic+1; i = ic
		dh  = yv(i,ny) - ycg(icG,nycg-1)
		Vf(ic,ny) = BCVal(3)*dh + vc(icG,nycg-1)
	END DO
 ELSEIF (BCType(3) .EQ. NOSLIP) THEN
	Vf(:,ny) = 0.D0
 END IF

 !---------------
 !Bottom Boundary
 !---------------
 IF (BCType(4) .EQ. INLET) THEN
	Vf(:,1) = BCVal(4)
 ELSEIF (BCType(4) .EQ. OUTLET) THEN
	DO ic = 1, nxc
		icG = ic+1; i = ic
		dh  = ycg(icG,2) - yv(i,1)
		Vf(ic,1) = BCVal(4)*dh + vc(icG,2)
	END DO
 ELSEIF (BCType(4) .EQ. NOSLIP) THEN
	Vf(:,1) = 0.D0
 END IF

 END SUBROUTINE SET_VBC_FC
 !----------------------------------------------------------------------------!


 SUBROUTINE SET_VBC_FC_SR

 USE GLOBAL
 USE MESH2D
 USE FLOW_PARAMS

 DOUBLE PRECISION :: dh, a, b, DBC, NBC

 INTEGER :: i, j, ic, jc, icG, jcG

 !--------------
 !Right Boundary
 !--------------
 IF (BCType(1) .EQ. INLET) THEN
	Ufsr(nx,:) = BCVal(1)
 ELSEIF (BCType(1) .EQ. OUTLET) THEN
	DO jc = 1, nyc
		jcG = jc+1; j = jc
		dh  = xv(nx,j) - xcg(nxcg-1,jcG)
		Ufsr(nx,jc) = BCVal(1)*dh + uc(nxcg-1,jcG)
	END DO
 ELSEIF (BCType(1) .EQ. NOSLIP) THEN
	Ufsr(nx,:) = 0.D0
 END IF

 !-------------
 !Left Boundary
 !-------------
 IF (BCType(2) .EQ. INLET) THEN
	Ufsr(1,:) = BCVal(2)
	DO jc = 1, nyc
		Uf(1,jc) = 1.D0 - 4.*BCVal(2)*yc(1,jc)**2
 	END DO
 ELSEIF (BCType(2) .EQ. OUTLET) THEN
	DO jc = 1, nyc
		jcG = jc+1; j = jc
		dh  = xcg(2,jcG) - xv(1,j)
		Ufsr(1,jc) = BCVal(2)*dh + uc(2,jcG)
	END DO
 ELSEIF (BCType(2) .EQ. NOSLIP) THEN
	Ufsr(1,:) = 0.D0
 END IF

 !------------
 !Top Boundary
 !------------
 IF (BCType(3) .EQ. INLET) THEN
	Vfsr(:,ny) = BCVal(3)
 ELSEIF (BCType(3) .EQ. OUTLET) THEN
	DO ic = 1, nxc
		icG = ic+1; i = ic
		dh  = yv(i,ny) - ycg(icG,nycg-1)
		Vfsr(ic,ny) = BCVal(3)*dh + vc(icG,nycg-1)
	END DO
 ELSEIF (BCType(3) .EQ. NOSLIP) THEN
	Vfsr(:,ny) = 0.D0
 END IF

 !---------------
 !Bottom Boundary
 !---------------
 IF (BCType(4) .EQ. INLET) THEN
	Vfsr(:,1) = BCVal(4)
 ELSEIF (BCType(4) .EQ. OUTLET) THEN
	DO ic = 1, nxc
		icG = ic+1; i = ic
		dh  = ycg(icG,2) - yv(i,1)
		Vfsr(ic,1) = BCVal(4)*dh + vc(icG,2)
	END DO
 ELSEIF (BCType(4) .EQ. NOSLIP) THEN
	Vfsr(:,1) = 0.D0
 END IF

 END SUBROUTINE SET_VBC_FC_SR
 !----------------------------------------------------------------------------!


 SUBROUTINE SET_VBC_VECTOR

 USE GLOBAL
 USE FLOW_PARAMS
 USE MESH2D
 IMPLICIT NONE

 DOUBLE PRECISION :: dxew, dyns
 DOUBLE PRECISION :: phiE, phiW, phiN, phiS
 DOUBLE PRECISION :: VBCD, VBCN, VAL
 DOUBLE PRECISION :: dh, a, b, DBC1, NBC1, DBC2, NBC2, test

 INTEGER :: i, j, ic, jc, icG, jcG, ICELL, JCELL
 INTEGER icE, icW, icN, icS

 VBC = 0.D0

 !Right boundary
 ic = nxc; icG = ic + 1

 IF (BCType(1) .EQ. INLET) THEN
	DO jc = 1, nyc
		jcG = jc + 1; j = jc
		ICELL = IC2CEL(ic,jc)

		dh	= xcg(nxcg,jcG)-xcg(nxcg-1,jcG)
		a	= xv(nx,j)-xcg(nxcg-1,jcG); b = xcg(nxcg,jcG) - xv(nx,j)

		DBC1= (BCVal(1)*dh - b*uc(nxcg-1,jcG))/a
		DBC2= -1.*vc(nxcg-1,jcG)

		VBC(ICELL,1) = VBC(ICELL,1) + 0.5D0*dt*Reinv*dxe2inv(ic,jc)*DBC1
		VBC(ICELL,2) = VBC(ICELL,2) + 0.5D0*dt*Reinv*dxe2inv(ic,jc)*DBC2
	END DO
 ELSEIF (BCType(1) .EQ. OUTLET) THEN
	DO jc = 1, nyc
		jcG = jc + 1; j = jc
		ICELL = IC2CEL(ic,jc)

		VBC(ICELL,1) = 0.D0!VBC(ICELL,1) - dxe2inv(ic,jc)*phiV(ICELL,1)
		VBC(ICELL,2) = 0.D0!VBC(ICELL,2) - dxe2inv(ic,jc)*phiV(ICELL,2)
		AmatV(ICELL,ICELL,1) = AmatV(ICELL,ICELL,1) - 0.5D0*dt*Reinv*dxe2inv(ic,jc)
		AmatV(ICELL,ICELL,2) = AmatV(ICELL,ICELL,2) - 0.5D0*dt*Reinv*dxe2inv(ic,jc)
	END DO
 ELSEIF (BCType(1) .EQ. NOSLIP) THEN
	DO jc = 1, nyc
		jcG = jc + 1; j = jc
		ICELL = IC2CEL(ic,jc)

		dh	= xcg(nxcg,jcG)-xcg(nxcg-1,jcG)
		a	= xv(nx,j)-xcg(nxcg-1,jcG); b = xcg(nxcg,jcG) - xv(nx,j)

		DBC1= -1.*uc(nxcg-1,jcG)
		DBC2= (BCVal(1)*dh - b*vc(nxcg-1,jcG))/a

		VBC(ICELL,1) = VBC(ICELL,1) + 0.5D0*dt*Reinv*dxe2inv(ic,jc)*DBC1
		VBC(ICELL,2) = VBC(ICELL,2) + 0.5D0*dt*Reinv*dxe2inv(ic,jc)*DBC2	
	END DO
 END IF

 !Left Boundary
 ic = 1; icG = ic + 1

 IF (BCType(2) .EQ. INLET) THEN
	DO jc = 1, nyc
		jcG = jc + 1; j = jc
		ICELL = IC2CEL(ic,jc)

		dh	= xcg(2,jcG)-xcg(1,jcG)
		a	= xcg(2,jcG)-xv(1,j); b = xv(1,j)-xcg(1,jcG)
		VAL = 1. - 4.*BCVal(2)*yc(ic,jc)**2

		DBC1= (VAL*dh - b*uc(2,jcG))/a
		DBC2= -1.*vc(2,jcG)

		VBC(ICELL,1) = VBC(ICELL,1) + 0.5D0*dt*Reinv*dxw2inv(ic,jc)*DBC1
		VBC(ICELL,2) = VBC(ICELL,2) + 0.5D0*dt*Reinv*dxw2inv(ic,jc)*DBC2
	END DO
 ELSEIF (BCType(2) .EQ. OUTLET) THEN
	DO jc = 1, nyc
		jcG = jc + 1; j = jc
		ICELL = IC2CEL(ic,jc)

!		VBC(ICELL,1) = VBC(ICELL,1) - dxw2inv(ic,jc)*phiV(ICELL,1)
!		VBC(ICELL,2) = VBC(ICELL,2) - dxw2inv(ic,jc)*phiV(ICELL,2)
		AmatV(ICELL,ICELL,1) = AmatV(ICELL,ICELL,1) - 0.5D0*dt*Reinv*dxw2inv(ic,jc)
		AmatV(ICELL,ICELL,2) = AmatV(ICELL,ICELL,2) - 0.5D0*dt*Reinv*dxw2inv(ic,jc)
	END DO
 ELSEIF (BCType(2) .EQ. NOSLIP) THEN
	DO jc = 1, nyc
		jcG = jc + 1; j = jc
		ICELL = IC2CEL(ic,jc)

		dh	= xcg(2,jcG)-xcg(1,jcG)
		a	= xcg(2,jcG)-xv(1,j); b = xv(1,j)-xcg(1,jcG)

		DBC1= -1.*uc(2,jcG)
		DBC2= (BCVal(2)*dh - b*vc(2,jcG))/a

		VBC(ICELL,1) = VBC(ICELL,1) + 0.5D0*dt*Reinv*dxw2inv(ic,jc)*DBC1
		VBC(ICELL,2) = VBC(ICELL,2) + 0.5D0*dt*Reinv*dxw2inv(ic,jc)*DBC2	
	END DO
 END IF

 !Top boundary
 jc = nyc; jcG = jc + 1

 IF (BCType(3) .EQ. INLET) THEN
	DO ic = 1, nxc
		icG = ic + 1; i = ic
		ICELL = IC2CEL(ic,jc)

		dh	= ycg(icG,nycg)-ycg(icG,nycg-1)
		a	= yv(i,ny)-ycg(icG,nycg-1); b = ycg(icG,nycg)-yv(i,ny)

		DBC1= -1.*uc(icG,nycg-1)
		DBC2= (BCVal(3)*dh - b*vc(icG,nycg-1))/a

		VBC(ICELL,1) = VBC(ICELL,1) + 0.5D0*dt*Reinv*dyn2inv(ic,jc)*DBC1
		VBC(ICELL,2) = VBC(ICELL,2) + 0.5D0*dt*Reinv*dyn2inv(ic,jc)*DBC2
	END DO
 ELSEIF (BCType(3) .EQ. OUTLET) THEN
	DO ic = 1, nxc
		icG = ic + 1
		ICELL = IC2CEL(ic,jc)

!		VBC(ICELL,1) = VBC(ICELL,1) - dyn2inv(ic,jc)*phiV(ICELL,1)
!		VBC(ICELL,2) = VBC(ICELL,2) - dyn2inv(ic,jc)*phiV(ICELL,2)
		AmatV(ICELL,ICELL,1) = AmatV(ICELL,ICELL,1) - 0.5D0*dt*Reinv*dyn2inv(ic,jc)
		AmatV(ICELL,ICELL,2) = AmatV(ICELL,ICELL,2) - 0.5D0*dt*Reinv*dyn2inv(ic,jc)
	END DO
 ELSEIF (BCType(3) .EQ. NOSLIP) THEN
	DO ic = 1, nxc
		icG = ic + 1; i = ic
		ICELL = IC2CEL(ic,jc)

		dh	= ycg(icG,nycg)-ycg(icG,nycg-1)
		a	= yv(i,ny)-ycg(icG,nycg-1); b = ycg(icG,nycg)-yv(i,ny)

		DBC1= (BCVal(3)*dh - b*uc(icG,nycg-1))/a
		DBC2= -1.*vc(icG,nycg-1)

		VBC(ICELL,1) = VBC(ICELL,1) + 0.5D0*dt*Reinv*dyn2inv(ic,jc)*DBC1
		VBC(ICELL,2) = VBC(ICELL,2) + 0.5D0*dt*Reinv*dyn2inv(ic,jc)*DBC2
	END DO
 END IF

 !Bottom boundary
 jc = 1; jcG = jc + 1

 IF (BCType(4) .EQ. INLET) THEN
	DO ic = 1, nxc
		icG = ic + 1; i = ic
		ICELL = IC2CEL(ic,jc)

		dh	= ycg(icG,2)-ycg(icG,1)
		a	= ycg(icG,2)-yv(i,1); b = yv(i,1)-ycg(icG,1) 

		DBC1= -1.*uc(icG,2)
		DBC2= (BCVal(4)*dh - b*vc(icG,2))/a

		VBC(ICELL,1) = VBC(ICELL,1) + 0.5D0*dt*Reinv*dys2inv(ic,jc)*DBC1
		VBC(ICELL,2) = VBC(ICELL,2) + 0.5D0*dt*Reinv*dys2inv(ic,jc)*DBC2
	END DO
 ELSEIF (BCType(3) .EQ. OUTLET) THEN
	DO ic = 1, nxc
		icG = ic + 1
		ICELL = IC2CEL(ic,jc)

!		VBC(ICELL,1) = VBC(ICELL,1) - dys2inv(ic,jc)*phiV(ICELL,1)
!		VBC(ICELL,2) = VBC(ICELL,2) - dys2inv(ic,jc)*phiV(ICELL,2)
		AmatV(ICELL,ICELL,1) = AmatV(ICELL,ICELL,1) - 0.5D0*dt*Reinv*dys2inv(ic,jc)
		AmatV(ICELL,ICELL,2) = AmatV(ICELL,ICELL,2) - 0.5D0*dt*Reinv*dys2inv(ic,jc)
	END DO
 ELSEIF (BCType(3) .EQ. NOSLIP) THEN
	DO ic = 1, nxc
		icG = ic + 1
		ICELL = IC2CEL(ic,jc)

		dh	= ycg(icG,2)-ycg(icG,1)
		a	= ycg(icG,2)-yv(i,1); b = yv(i,1)-ycg(icG,1) 

		DBC1= (BCVal(4)*dh - b*vc(icG,2))/a
		DBC2= -1.*uc(icG,2)

		VBC(ICELL,1) = VBC(ICELL,1) + 0.5D0*dt*Reinv*dys2inv(ic,jc)*DBC1
		VBC(ICELL,2) = VBC(ICELL,2) + 0.5D0*dt*Reinv*dys2inv(ic,jc)*DBC2
	END DO
 END IF

! WRITE(*,'(A)') ' DONE!'
! PRINT*,''

 END SUBROUTINE SET_VBC_VECTOR
 !----------------------------------------------------------------------------!
