
 SUBROUTINE ALLOCATE_MEMORY

 USE GLOBAL, ONLY: iIntlBody
 IMPLICIT NONE

 WRITE(*,'(A)',ADVANCE='NO') ' ALLOCATING GRID MEMORY ...'
 CALL ALLOCATE_GRID_MEMORY
 WRITE(*,'(A)') ' DONE!'
 WRITE(*,'(A)') ''

 WRITE(*,'(A)',ADVANCE='NO') ' ALLOCATING SOLUTION VARIABLE MEMORY ...'
 CALL ALLOCATE_FLOW_MEMORY
 WRITE(*,'(A)') ' DONE!'
 PRINT*,''

 WRITE(*,'(A)',ADVANCE='NO') ' ALLOCATING POISSON MEMORY ...'
 CALL ALLOCATE_POISSON_MEMORY
 WRITE(*,'(A)') ' DONE!'
 PRINT*,''

 IF (iIntlBody) THEN
	 WRITE(*,'(A)',ADVANCE='NO') ' ALLOCATING IBM MEMORY ...'
	 CALL ALLOCATE_IBM_MEMORY
	 WRITE(*,'(A)') ' DONE!'
	 PRINT*,''
 END IF

 END SUBROUTINE ALLOCATE_MEMORY
 !----------------------------------------------------------------------------!


 SUBROUTINE ALLOCATE_GRID_MEMORY

 USE GLOBAL
 USE MESH2D
 IMPLICIT NONE

 ALLOCATE(xv(nx,ny), yv(nx,ny))
 ALLOCATE(xvg(nxg,nyg),yvg(nxg,nyg))
 ALLOCATE(xc(nxc,nyc),yc(nxc,nyc))
 ALLOCATE(xcg(nxcg,nycg),ycg(nxcg,nycg))

 ALLOCATE(dxcinv(nxc,nyc),dycinv(nxc,nyc))
 ALLOCATE(dxcinv_NC(nxc,nyc),dycinv_NC(nxc,nyc))
 ALLOCATE(dxfinv(nx,nyc),dyfinv(nxc,ny))
 ALLOCATE(dxe2inv(nxc,nyc),dxw2inv(nxc,nyc),dyn2inv(nxc,nyc),dys2inv(nxc,nyc))

 ALLOCATE(IC2CEL(nxc,nyc),ICEL2C(NCELL,2),IC2CELG(nxcg,nycg),ICELG2C(NCELLG,2))
 ALLOCATE(IC2NB(nxc,nyc,4),IC2NBG(nxcg,nycg,4))
 ALLOCATE(ibE(nxc,nyc),ibW(nxc,nyc),ibN(nxc,nyc),ibS(nxc,nyc))

 END SUBROUTINE ALLOCATE_GRID_MEMORY
 !----------------------------------------------------------------------------!


 SUBROUTINE ALLOCATE_FLOW_MEMORY

 USE GLOBAL
 USE MESH2D
 USE FLOW_PARAMS
 IMPLICIT NONE

 ALLOCATE(Uvar(nxc,nyc),f_Uvar(nxc,nyc))
 
 !For NS-Solver
 ALLOCATE(uc(nxcg,nycg), vc(nxcg,nycg), pc(nxcg,nycg))
 ALLOCATE(usr(nxcg,nycg), vsr(nxcg,nycg))
 ALLOCATE(Uf(nx,nyc), Vf(nxc,ny))
 ALLOCATE(Ufsr(nx,nyc), Vfsr(nxc,ny))

 ALLOCATE(FconU(nxcg,nycg),FconV(nxcg,nycg))
 ALLOCATE(FconU_old(nxcg,nycg),FconV_old(nxcg,nycg))
 ALLOCATE(FdifU(nxcg,nycg),FdifV(nxcg,nycg))

! ALLOCATE(AmatV(NCELL,NCELL,2))
 ALLOCATE(rhsV(NCELL,2),phiV(NCELL,2),VBC(NCELL,2))
 ALLOCATE(AmatVG(NCELLG,NCELLG),AmatPG(NCELLG,NCELLG))

 END SUBROUTINE ALLOCATE_FLOW_MEMORY
 !----------------------------------------------------------------------------!


 SUBROUTINE ALLOCATE_IBM_MEMORY

 USE GLOBAL
 USE MESH2D
 USE IBM_PARAMS
 IMPLICIT NONE

 !For Internal Body flow simulations
 ALLOCATE(IBLANK(nxcg,nycg),GMAP(nxcg,nycg))

 END SUBROUTINE ALLOCATE_IBM_MEMORY
 !----------------------------------------------------------------------------!


 SUBROUTINE FREE_MEMORY

 USE GLOBAL 
 USE MESH2D
 USE POISSON_PARAMS
 USE FLOW_PARAMS
 USE IBM_PARAMS
 IMPLICIT NONE

 WRITE(*,'(A)',ADVANCE='NO') ' FREEING MEMORY ...'

 !-----------------------------------------------!
 !Deallocate Grid params:
 !----------------------
 DEALLOCATE(xv,yv,xc,yc)
 DEALLOCATE(xvg,yvg,xcg,ycg)
 DEALLOCATE(dxcinv,dycinv)
 DEALLOCATE(dxcinv_NC,dycinv_NC)
 DEALLOCATE(dxfinv,dyfinv)
 DEALLOCATE(dxe2inv,dxw2inv,dyn2inv,dys2inv)
 DEALLOCATE(IC2CEL,ICEL2C,IC2NB,IC2NBG)
 DEALLOCATE(ibE,ibW,ibN,ibS)
 !-----------------------------------------------!

 !-----------------------------------------------!
 !Deallocate Poisson params:
 !-------------------------
 DEALLOCATE(phi,Uvar,f_Uvar)
 !-----------------------------------------------!

 !-----------------------------------------------!
 !Deallocate flow vars:
 !--------------------
 DEALLOCATE(uc,vc,pc)
 DEALLOCATE(Uf,Vf)
 DEALLOCATE(usr,vsr)
 DEALLOCATE(Ufsr,Vfsr)
 DEALLOCATE(FconU,FconV)
 DEALLOCATE(FconU_old,FconV_old)
 DEALLOCATE(FdifU,FdifV)
! DEALLOCATE(AmatV)
 DEALLOCATE(rhsV,phiV,VBC)
 DEALLOCATE(AmatVG,AmatPG)
 DEALLOCATE(xProbe,yProbe,uProbe,vProbe,pProbe)
 DEALLOCATE(LOC_PROBE)
 !-----------------------------------------------!

 !-----------------------------------------------!
 !Deallocate IBM vars:
 !--------------------
 IF (iIntlBody) THEN
	DEALLOCATE(IBLANK,GMAP,IS2C)
	DEALLOCATE(IGC2C,IGC2FNB)
 END IF

 WRITE(*,'(A)') ' DONE!'
 WRITE(*,'(A)') ''
 
 END SUBROUTINE FREE_MEMORY
 !----------------------------------------------------------------------------!
