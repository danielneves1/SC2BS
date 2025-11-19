
!**************************************************************************/
Subroutine time_verb(ds)
	
!**************************************************************************/

	!Shows in either a file (232) (ds==1) or in terminal (*) (ds!=1) the expected run time and its current status

use variables
implicit none

integer, intent(in) :: ds

dum = sum(timetest)/it *(lent-it-1)



if (ds==1) then

	if( dum < 60) then
	    
		write(232, "(A, F3.0, A, I0,A, I0,A)", advance="no") char(13) // "Approx Time left = ", dum ," secs | Iteration: ", it, "/", lent-1

	endif

	if( dum<3600 .and. dum>60) then
		dummy = floor( dum/60 )
		dummy2 = ( dum/60-dummy)*60


		write(232, "(A, I0, A,F3.0, A,I0, A,I0 ,A)", advance="no") char(13) // "Approx Time left = ",  dummy ,":", dummy2 ," min | Iteration: ", it, "/", lent-1

	endif

	if(dum>3600) then

		dummy = floor( dum/3600 )
		dummy3 = floor( (dum/3600-dummy)*60 )
		dummy2 = (  ( dum/3600-dummy )*60 -dummy3  )*60
		write(*, "(A, I0, A,I0, A,F3.0, A, I0,A, I0,A)", advance="no") char(13) // "Approx Time left  = ", dummy ,":", dummy3,":",dummy2," hour | Iteration: ", it, "/", lent-1

	endif

else

	if( dum < 60) then
	    
		write(*, "(A, F3.0, A, I0,A, I0,A)", advance="no") char(13) // "Approx Time left = ", dum ," secs | Iteration: ", it, "/", lent-1

	endif

	if( dum<3600 .and. dum>60) then
		dummy = floor( dum/60 )
		dummy2 = ( dum/60-dummy)*60


		write(*, "(A, I0, A,F3.0, A,I0, A,I0 ,A)", advance="no") char(13) // "Approx Time left = ",  dummy ,":", dummy2 ," min | Iteration: ", it, "/", lent-1

	endif

	if(dum>3600) then

		dummy = floor( dum/3600 )
		dummy3 = floor( (dum/3600-dummy)*60 )
		dummy2 = (  ( dum/3600-dummy )*60 -dummy3  )*60
		write(*, "(A, I0, A,I0, A,F3.0, A, I0,A, I0,A)", advance="no") char(13) // "Approx Time left  = ", dummy ,":", dummy3,":",dummy2," hour | Iteration: ", it, "/", lent-1

	endif



endif

end subroutine time_verb



!**************************************************************************/
Subroutine writer(full)

		! full == 0 --> writes only the instants where alpha ratio are frac apart (can be changed)
		!	      These instant are WRITTEN in (66)	(as "tout.txt"), not to be mistaken by 
		!	      (6666) (as "instants.txt") which are ALL the EVALUATED instants.	

		! full != 0 --> writes every instant (frac will not matter). And "instants.txt" should be equiv. to "tout.txt".
!**************************************************************************/
use variables
implicit none

integer :: full

if (full == 0) then 

!Only writes ratios that are frac apart

	if (alfaratio_t(it)<1.01*(alfaratio_t(1)-frac*n) .and. alfaratio_t(it)>0.99*(alfaratio_t(1)-frac*n)) then
		do i=1,lenx
			write(3,*) psi_k(i)
			write(2,*) V_k(i)
		enddo
		write(66,*) timearray(it)
		n=n+1
	endif

	if (timearray(it)-to-tau<0) then
		m=9
	else
		m=m+1
	endif

!After the external potential evaporates prints all iterations

	if (timearray(it)-to-tau>=0 .and. m==10) then
		m=0
		do i=1,lenx
			write(3,*) psi_k(i)
			write(2,*) V_k(i)
		enddo
		write(66,*) timearray(it)
	endif

else

!Write all iterations

	do i=1,lenx
		write(3,*) psi_k(i)
		write(2,*) V_k(i)
	enddo
		write(66,*) timearray(it)
	

endif


end subroutine writer