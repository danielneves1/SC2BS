!*************************************************************************/
  SUBROUTINE CloudPotential(V, psi)
!**************************************************************************/


Use Variables
IMPLICIT NONE
Complex*16, dimension(lenx),intent(in) :: psi
REAL*8, dimension(lenx), intent(out):: V



!V(lenx) = 1

!V(lenx-1) = 1 - dxx / x(lenx)

!do i=lenx-1, 2, -1

!V(i-1) =  (-V(i)*diag(i) - V(i+1)*topdiag(i) - abs(psi(i))**2/x(i) ) /subdiag(i)
!enddo


!***********************************************************************************************************
 !solves Laplacian(V_k) = -|psi|^2/r for a Laplacian written as tridiagonal matrix as in simple finite diff.
!***********************************************************************************************************

call solve_tridiagr(subdiag, diag, topdiag, -abs(psi)**2/x, V, dxx(lenx) / x(lenx),  lenx)

End Subroutine CloudPotential
