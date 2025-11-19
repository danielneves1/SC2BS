



!*************************************************************************/
  SUBROUTINE tempos(c, t, kk)
!**************************************************************************/

! Subroutine to obtain the simulation instants based on Hawking evaporation time scale
! kk is the input which makes alpha to only vary by kk amount
 
use Variables
Implicit none
	real*8 ,intent(in) :: t, kk
	real*8  :: r
	real*8 , intent(out)  :: c
	r=(t-to)/tau
	
	if(r<1) then
		c=tau*(1 - ( (1-r)**(1d0/3d0)-kk )**(3d0))
		if(c-(t-to)>0.05) then
			c=t-to+0.05
		endif
	else
		c=1E38
	endif
	
	

End Subroutine tempos


