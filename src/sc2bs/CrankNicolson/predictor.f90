!*************************************************************************/
  SUBROUTINE predictor_Article()
!**************************************************************************/

use Variables
Implicit none

call wf(psi_k1, V_k/x + pot_r)
Call CloudPotential(V_k1, psi_k1)

V0 = (V_k1 + V_k)/2
pot_r = (alfaratio_t(it)+alfaratio_t(it+1))/2/x

call wf(psi_k2, V0/x + pot_r)
call CloudPotential(V_k2, psi_k2)



End Subroutine predictor_Article