 SUBROUTINE ITERATE_JACOBI(iter_hold)

 USE POISSON_PARAMS
 USE MESH2D
 IMPLICIT NONE

 DOUBLE PRECISION, ALLOCATABLE :: DI(:), LU(:,:), LU_phi(:), phi_old(:)

 INTEGER, INTENT(IN) :: iter_hold
 INTEGER :: i, j, iter_p, ICELL

 !Allocate memory
 ALLOCATE(DI(NCELL))		!D^-1
 ALLOCATE(LU(NCELL,NCELL))	!(L+U)
 ALLOCATE(LU_phi(NCELL))	!(L+U)*phi
 ALLOCATE(phi_old(NCELL))	

 iter_p = 0; pflag = .FALSE.
 LU = 0.D0; LU_phi = 0.D0; DI = 0.D0

 !Arrange coefficient matrices
 DO ICELL = 1, NCELL
	DI(ICELL)   = 1.D0/Amat(ICELL,ICELL)
	LU(ICELL,:) = Amat(ICELL,:)
	LU(ICELL,ICELL) = 0.D0
 END DO

 !Iter_hold provides the solver to run up to a prescribed iteration count
 !Iter_hold is used to run the smoother for Multigrid
 DO WHILE ( (iter_p .LT. iter_hold .AND. (.NOT. pflag)) &
	   .OR. (iter_p .LT. PITERMAX  .AND. (.NOT. pflag)) )

	phi_old = phi

	CALL SET_PBC

	CALL MATPROD(LU,phi_old,LU_phi,NCELL)

	!Compute new phi
	phi = DI*(rhs + PBC - LU_phi)

	WRITE(*,'(A,I5)') ' >> ITER: ', iter_p

	CALL COMP_ERROR(iter_p)

	iter_p = iter_p + 1

 END DO

 !Free memory
 DEALLOCATE(DI,LU,LU_phi,phi_old)

 END SUBROUTINE ITERATE_JACOBI
 !----------------------------------------------------------------------------!


 SUBROUTINE ITERATE_SOR(iter_hold)

 USE MESH2D
 USE POISSON_PARAMS
 IMPLICIT NONE

 DOUBLE PRECISION, ALLOCATABLE :: temp(:), DwL(:,:), wUw1D(:,:), wUw1D_phi(:), phi_old(:)

 INTEGER, INTENT(IN) :: iter_hold
 INTEGER :: i, j, iter_p, ICELL

 iter_p = 1; pflag = .FALSE.

 !Allocate memory
 ALLOCATE(DwL(NCELL,NCELL))		!(D+w*L)
 ALLOCATE(wUw1D(NCELL,NCELL))	!(wU+(w-1)*D)
 ALLOCATE(wUw1D_phi(NCELL))		!(wU+(w-1)*D)phi
 ALLOCATE(phi_old(NCELL))
 ALLOCATE(temp(NCELL))

 !Arrange coefficient matrices
 DO ICELL = 1, NCELL
	DwL(ICELL,1:ICELL-1)	= omega*Amat(ICELL,1:ICELL-1)
	DwL(ICELL,ICELL)		= Amat(ICELL,ICELL)

	wUw1D(ICELL,ICELL+1:)	= omega*(Amat(ICELL,ICELL+1:))
	wUw1D(ICELL,ICELL)		= (omega-1.D0)*Amat(ICELL,ICELL)
 END DO

 !Iter_hold provides the solver to run up to a prescribed iteration count
 !Iter_hold is used to run the smoother for Multigrid
 DO WHILE ( (iter_p .LT. iter_hold .AND. (.NOT. pflag)) )

	phi_old = phi

	CALL SET_PBC

	CALL MATPROD(wUw1D,phi_old,wUw1D_phi,NCELL)

	temp = omega*(rhs + PBC) - wUw1D_phi

	!Compute new phi
	CALL INV_FWSUB_SOLVE(DwL,temp,phi,NCELL)

!	WRITE(*,'(A,I5)') ' >> ITER: ', iter_p

	CALL COMP_ERROR(iter_p)

	iter_p = iter_p + 1

 END DO

 !Free memory
 DEALLOCATE(DwL,wUw1D,wUw1D_phi,phi_old,temp)

 END SUBROUTINE ITERATE_SOR
 !----------------------------------------------------------------------------!
 

 SUBROUTINE ITERATE_CG(iter_hold)

 USE MESH2D
 USE POISSON_PARAMS
 IMPLICIT NONE

 INTERFACE
	DOUBLE PRECISION FUNCTION DOTPROD(v1,v2,LENGTH)
	
	DOUBLE PRECISION, INTENT(IN) :: v1(LENGTH), v2(LENGTH)

	INTEGER, INTENT(IN) :: LENGTH

	END FUNCTION DOTPROD
 END INTERFACE

 DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: Amat2
 DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: PBC2, DI
 DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: b, r, rnew, rtil, p, Ap
 DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: phi_old, Aphi0
 DOUBLE PRECISION :: alpha, beta

 INTEGER, INTENT(IN) :: iter_hold
 INTEGER :: ic, jc, ICELL
 INTEGER :: icE, icW, icN, icS
 INTEGER :: iter_p, i

 ALLOCATE(Amat2(NCELL,NCELL), PBC2(NCELL))
 ALLOCATE(b(NCELL))
 ALLOCATE(r(NCELL),rnew(NCELL),rtil(NCELL))
 ALLOCATE(p(NCELL),Ap(NCELL))
 ALLOCATE(phi_old(NCELL),Aphi0(NCELL))

 !Set boundary condition
 CALL SET_PBC

 !--------------------------------------------------------------!
 !Note:
 !In this implementation, the boundary condition vector (PBC) is 
 !reflected only in the initial computation. So the A matrix is 
 !modified to include only the Neumann BC. This component is sub-
 !tracted from the PBC vector.
 !--------------------------------------------------------------!

 !Set modified A matrix and PBC vector
 Amat2 = Amat; PBC2 = PBC
 DO ICELL = 1, NCELL
	ic = ICEL2C(ICELL,1); jc = ICEL2C(ICELL,2)

	icE = IC2NB(ic,jc,1)
	icW = IC2NB(ic,jc,2)
	icN = IC2NB(ic,jc,3)
	icS = IC2NB(ic,jc,4)

	!Add Neumann BC to A matrix
	Amat2(ICELL,ICELL) = Amat2(ICELL,ICELL) + &
			  ( REAL(ibE(ic,jc)*(1-PBCType(1)),KIND=8)*dxe2inv(ic,jc) + &
				REAL(ibW(ic,jc)*(1-PBCType(2)),KIND=8)*dxw2inv(ic,jc) + &
				REAL(ibN(ic,jc)*(1-PBCType(3)),KIND=8)*dyn2inv(ic,jc) + &
				REAL(ibS(ic,jc)*(1-PBCType(4)),KIND=8)*dys2inv(ic,jc) )

	!Remove Neumann BC from PBC vector
	PBC2(ICELL) = PBC2(ICELL) +	&
			  ( REAL(ibE(ic,jc)*(1-PBCType(1)),KIND=8)*dxe2inv(ic,jc) + &
				REAL(ibW(ic,jc)*(1-PBCType(2)),KIND=8)*dxw2inv(ic,jc) + &
				REAL(ibN(ic,jc)*(1-PBCType(3)),KIND=8)*dyn2inv(ic,jc) + &
				REAL(ibS(ic,jc)*(1-PBCType(4)),KIND=8)*dys2inv(ic,jc) )*(phi(ICELL))
 END DO

 !Initialize computation vectors.
 CALL MATPROD(Amat2,phi,Aphi0,NCELL)
 rnew	= rhs + PBC2 - Aphi0
 p		= rnew

 iter_p = 0; pflag = .FALSE.

 !Iter_hold provides the solver to run up to a prescribed iteration count
 !Iter_hold is used to run the smoother for Multigrid
 DO WHILE ( (iter_p .LT. iter_hold .AND. (.NOT. pflag)) )

	phi_old = phi; r = rnew

	CALL MATPROD(Amat2,p,Ap,NCELL)

	alpha = DOTPROD(r,r,NCELL)/DOTPROD(p,Ap,NCELL)

	phi  = phi_old	+ alpha*p
	rnew = r 		- alpha*Ap

	beta = DOTPROD(rnew,rnew,NCELL)/DOTPROD(r,r,NCELL)
	p 	 = rnew + beta*p

!	WRITE(*,'(A,I5)') ' >> ITER: ', iter_p

	CALL COMP_ERROR(iter_p)

	iter_p = iter_p + 1
 END DO

 DEALLOCATE(Amat2,PBC2)
 DEALLOCATE(b,r,rnew,p,Ap,phi_old) 

 END SUBROUTINE ITERATE_CG
 !----------------------------------------------------------------------------!


 SUBROUTINE INV_FWSUB_SOLVE (L,b,x,LENGTH)

 IMPLICIT NONE

 DOUBLE PRECISION, INTENT(IN)  :: L(LENGTH,LENGTH), b(LENGTH)
 DOUBLE PRECISION, INTENT(OUT) :: x(LENGTH)
 DOUBLE PRECISION :: LI(LENGTH,LENGTH), temp(LENGTH), summ

 INTEGER :: i, j, m, LENGTH
 
 !Note that L in this case includes the diagonal entries
 temp = 0.D0

 temp(1) = b(1)/L(1,1)

 DO m = 2, LENGTH
	DO i = 1, m-1
		temp(m) = temp(m) + L(m,i)*temp(i)
	END DO

	temp(m) = (b(m)-temp(m))/L(m,m)
 END DO

 x = temp

 END SUBROUTINE INV_FWSUB_SOLVE
 !----------------------------------------------------------------------------!
 

 SUBROUTINE MATPROD(A,b,Ab,LENGTH)

 IMPLICIT NONE

 DOUBLE PRECISION, INTENT(IN) :: A(LENGTH,LENGTH), b(LENGTH)
 DOUBLE PRECISION :: Ab(LENGTH), summ

 INTEGER :: i, j, LENGTH

 Ab = 0.0

 DO i = 1, LENGTH
	DO j = 1, LENGTH
		Ab(i) = Ab(i) + A(i,j)*b(j)
	END DO
 END DO 

 END SUBROUTINE MATPROD
 !----------------------------------------------------------------------------!
 

 DOUBLE PRECISION FUNCTION DOTPROD(v1, v2, LENGTH)

 DOUBLE PRECISION, INTENT(IN) :: v1(LENGTH), v2(LENGTH)
 DOUBLE PRECISION :: summ

 INTEGER, INTENT(IN) :: LENGTH
 INTEGER :: i

 summ = 0.D0

 DO i = 1, LENGTH
	summ = summ + v1(i)*v2(i)
 END DO

 DOTPROD = summ

 END FUNCTION DOTPROD
 !----------------------------------------------------------------------------!
