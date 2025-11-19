
!*************************************************************************/
  SUBROUTINE WF(psi, pot_term)
!*************************************************************************/

!# solves  A.u_k1 = B.u_k for u_k1 where A and B are tridiagonal matrices

Use Variables
IMPLICIT NONE

real*8, dimension(lenx), intent(in) :: pot_term
Complex*16, dimension(lenx), intent(out) :: psi

!# compute B.u_k
Call tridiag_mult( -subdiag/2, timediag - diag/2 - pot_term/2 , -topdiag/2, lenx, psi_k, Bpsi)

!# computes u_k,1 through A inversion
call solve_tridiag(subdiag/2, timediag + diag/2 + pot_term/2, topdiag/2, Bpsi, psi, lenx)


End Subroutine wf



