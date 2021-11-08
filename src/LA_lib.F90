 SUBROUTINE ITERATE_JACOBI(Amat, phi, rhs, iter_hold, L2NORM)

 USE GLOBAL
 USE MESH2D
 IMPLICIT NONE

 DOUBLE PRECISION, INTENT(IN) :: Amat(NCELL,NCELL), rhs(NCELL)
 DOUBLE PRECISION, INTENT(INOUT) :: phi(NCELL)
 DOUBLE PRECISION, INTENT(OUT) :: L2NORM
 DOUBLE PRECISION, ALLOCATABLE :: DI(:), LU(:,:), LU_phi(:), phi_old(:)

 INTEGER, INTENT(IN) :: iter_hold
 INTEGER :: i, j, iter_p, ICELL

 LOGICAL :: flag

 !Allocate memory
 ALLOCATE(DI(NCELL))		!D^-1
 ALLOCATE(LU(NCELL,NCELL))	!(L+U)
 ALLOCATE(LU_phi(NCELL))	!(L+U)*phi
 ALLOCATE(phi_old(NCELL))	

 iter_p = 0; flag = .FALSE.
 LU = 0.D0; LU_phi = 0.D0; DI = 0.D0

 !Arrange coefficient matrices
 DO ICELL = 1, NCELL
	DI(ICELL)   = 1.D0/Amat(ICELL,ICELL)
	LU(ICELL,:) = Amat(ICELL,:)
	LU(ICELL,ICELL) = 0.D0
 END DO

 !Iter_hold provides the solver to run up to a prescribed iteration count
 !Iter_hold is used to run the smoother for Multigrid
 DO WHILE (iter_p .LT. iter_hold)

	phi_old = phi

	CALL MATPROD(LU,phi_old,LU_phi,NCELL)

	!Compute new phi
	phi = DI*(rhs - LU_phi)

!	WRITE(*,'(A,I5)') ' >> ITER: ', iter_p

	CALL COMP_ERROR(Amat, phi, rhs, iter_p, L2NORM)

	iter_p = iter_p + 1

 END DO

 !Free memory
 DEALLOCATE(DI,LU,LU_phi,phi_old)

 END SUBROUTINE ITERATE_JACOBI
 !----------------------------------------------------------------------------!


 SUBROUTINE ITERATE_SOR(Amat, phi, rhs, iter_hold, L2NORM)

 USE GLOBAL
 USE MESH2D
 IMPLICIT NONE

 DOUBLE PRECISION, INTENT(IN) :: Amat(NCELL,NCELL), rhs(NCELL)
 DOUBLE PRECISION, INTENT(INOUT) :: phi(NCELL)
 DOUBLE PRECISION, INTENT(OUT) :: L2NORM
 DOUBLE PRECISION, ALLOCATABLE :: temp(:), DwL(:,:), wUw1D(:,:), wUw1D_phi(:), phi_old(:)

 INTEGER, INTENT(IN) :: iter_hold
 INTEGER :: i, j, iter_p, ICELL

 LOGICAL :: flag

 iter_p = 1; flag = .FALSE.

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
 DO WHILE (iter_p .LE. iter_hold)

	phi_old = phi

	CALL MATPROD(wUw1D,phi_old,wUw1D_phi,NCELL)

	temp = omega*rhs - wUw1D_phi

	!Compute new phi
	CALL INV_FWSUB_SOLVE(DwL,temp,phi,NCELL)

!	WRITE(*,'(A,I5)') ' >> ITER: ', iter_p

	CALL COMP_ERROR(Amat, phi, rhs, iter_p, L2NORM)

	iter_p = iter_p + 1

 END DO

 !Free memory
 DEALLOCATE(DwL,wUw1D,wUw1D_phi,phi_old,temp)

 END SUBROUTINE ITERATE_SOR
 !----------------------------------------------------------------------------!
 

 SUBROUTINE ITERATE_CG(Amat, phi, rhs, iter_hold, L2NORM)

 USE MESH2D
 USE GLOBAL
 IMPLICIT NONE

 INTERFACE
	DOUBLE PRECISION FUNCTION DOTPROD(v1,v2,LENGTH)
	
	DOUBLE PRECISION, INTENT(IN) :: v1(LENGTH), v2(LENGTH)

	INTEGER, INTENT(IN) :: LENGTH

	END FUNCTION DOTPROD
 END INTERFACE

 DOUBLE PRECISION, INTENT(IN) :: Amat(NCELL,NCELL), rhs(NCELL)
 DOUBLE PRECISION, INTENT(INOUT) :: phi(NCELL)
 DOUBLE PRECISION, INTENT(OUT) :: L2NORM
 DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: Amat2
 DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: PBC2, DI
 DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: b, r, rnew, rtil, p, Ap
 DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: phi_old, Aphi0
 DOUBLE PRECISION :: alpha, beta

 INTEGER, INTENT(IN) :: iter_hold
 INTEGER :: ic, jc, ICELL
 INTEGER :: icE, icW, icN, icS
 INTEGER :: iter_p, i

 LOGICAL :: flag

 ALLOCATE(Amat2(NCELL,NCELL), PBC2(NCELL))
 ALLOCATE(b(NCELL))
 ALLOCATE(r(NCELL),rnew(NCELL),rtil(NCELL))
 ALLOCATE(p(NCELL),Ap(NCELL))
 ALLOCATE(phi_old(NCELL),Aphi0(NCELL))

 !Initialize computation vectors.
 CALL MATPROD(Amat,phi,Aphi0,NCELL)
 rnew	= rhs - Aphi0
 p		= rnew

 iter_p = 0; flag = .FALSE.

 !Iter_hold provides the solver to run up to a prescribed iteration count
 !Iter_hold is used to run the smoother for Multigrid
 DO WHILE (iter_p .LE. iter_hold)

	phi_old = phi; r = rnew

	CALL MATPROD(Amat,p,Ap,NCELL)

	alpha = DOTPROD(r,r,NCELL)/DOTPROD(p,Ap,NCELL)

	phi  = phi_old	+ alpha*p
	rnew = r 		- alpha*Ap

	beta = DOTPROD(rnew,rnew,NCELL)/DOTPROD(r,r,NCELL)
	p 	 = rnew + beta*p

	CALL COMP_ERROR(Amat, phi, rhs, iter_p, L2NORM)

	iter_p = iter_p + 1
 END DO

 DEALLOCATE(Amat2,PBC2)
 DEALLOCATE(b,r,rnew,p,Ap,phi_old) 

 END SUBROUTINE ITERATE_CG
 !----------------------------------------------------------------------------!
 

 SUBROUTINE COMP_ERROR(Amat, phi, rhs, iter_p, L2NORM)

 USE MESH2D
 USE POISSON_PARAMS, ONLY: P_tol
 IMPLICIT NONE
 
 DOUBLE PRECISION, INTENT(IN) :: Amat(NCELL,NCELL), rhs(NCELL)
 DOUBLE PRECISION, INTENT(INOUT) :: phi(NCELL)
 DOUBLE PRECISION, INTENT(OUT) :: L2NORM
 DOUBLE PRECISION :: ERROR, ESUM
 DOUBLE PRECISION, ALLOCATABLE :: A_phi(:)

 INTEGER, INTENT(IN) :: iter_p
 INTEGER :: ic, jc, ICELL, counter

 LOGICAL :: flag

 ALLOCATE(A_phi(NCELL))

 ESUM = 0.D0;

 CALL MATPROD(Amat,phi,A_phi,NCELL)

 !Assemble coefficient and rhs for interior points
 DO ICELL = 1, NCELL
	ic = ICEL2C(ICELL,1); jc = ICEL2C(ICELL,2)
	ERROR = rhs(ICELL) - A_phi(ICELL)
	ESUM  = ESUM + ERROR**2

	counter = counter + 1
 END DO

 L2NORM = SQRT(ESUM/counter)

!	WRITE(*,'(A,E15.8)') '    L2NORM: ', L2NORM
!	WRITE(*,'(A)') ''

! IF (L2NORM .LT. P_tol) THEN
!	flag = .TRUE. 
!	WRITE(*,'(A)') ''
!	WRITE(*,'(A,I4)')    '    SOLUTION CONVERGED AFTER ITERATION: ', iter_p
!	WRITE(*,'(A,E15.8)') '    L2NORM: ', L2NORM
! END IF

 DEALLOCATE(A_phi)

 END SUBROUTINE COMP_ERROR
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


 SUBROUTINE TDMA(Mat,d,res,LENGTH)

 IMPLICIT NONE

 DOUBLE PRECISION :: Mat(LENGTH,LENGTH)
 DOUBLE PRECISION, INTENT(IN) :: d(LENGTH)
 DOUBLE PRECISION, INTENT(OUT) :: res(LENGTH)
 DOUBLE PRECISION, DIMENSION(LENGTH) :: a, b, c, cp, dp, temp

 INTEGER, INTENT(IN) :: LENGTH
 INTEGER :: n, i, j

 a = 0.; b = 0.; c =0.

 !Assign a, b, c
 DO i = 1, LENGTH
	IF (i .NE. 1)		a(i) = Mat(i,i-1)
						b(i) = Mat(i,i)
	IF (i .NE. LENGTH)	c(i) = Mat(i,i+1)
 END DO

 !Forward sweep
 DO i = 1, LENGTH
	IF (i .EQ. 1) THEN
		cp(i) = c(i)/b(i)
		dp(i) = d(i)/b(i)
	ELSEIF (i .EQ. LENGTH) THEN
		dp(i) = ( d(i)-a(i)*dp(i-1) )/( b(i)-a(i)*cp(i-1) )
	ELSE
		cp(i) = c(i)/( b(i)-a(i)*cp(i-1) )
		dp(i) = ( d(i)-a(i)*dp(i-1) )/( b(i)-a(i)*cp(i-1) )
	END IF
 END DO

 !Back substitution
 temp(LENGTH) = dp(LENGTH)
 DO i = LENGTH-1,1,-1
	temp(i) = dp(i) - cp(i)*temp(i+1)
 END DO
 res = temp

 END SUBROUTINE TDMA
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
