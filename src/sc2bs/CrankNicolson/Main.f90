
!***************************************************
PROGRAM CrankNicolson
!***************************************************


USE Variables

IMPLICIT NONE

!**********************************************************************************************
!**************************** Set up the output and input file names* *************************
!**********************************************************************************************

Ficheiro = "Input/inp.txt"
potout = "Output/Potentialout.txt"
wfin = "Input/WF.txt"
wfout = "Output/wfout.txt"
xout="Output/xout.txt"
xpoints ="Input/xpoints.txt"
tout ="Output/tout.txt"
Energ ="Output/energy.txt"
Norma ="Output/norm.txt"

open(4321,file="last_wf.txt")


OPen(232, File="TimeLeft.txt")
OPEN(1, FILE=Ficheiro)
Open(2, File=potout)
Open(3, File=wfout)
Open(4, File=xpoints)
Open(66, File=tout)
Open(8, File=energ)
Open(9, File=norma)

open(22, file= xout)
Open(70, File=wfin)
open(6666, file="Output/instants.txt")

!**********************************************************************************************
!**********************************************************************************************

!**********************************************************************************************
!**************************** initialize the simulation parameters ****************************
!**********************************************************************************************

	call Initialize
!**********************************************************************************************
!**********************************************************************************************


!**********************************************************************************************
!**********************************************************************************************
	!Compute the initial energy 

call CloudPotential(V_k, psi0)

call gradc(wf_grad, psi0/x, x, lenx )

call integrate((abs(wf_grad)**2-(V_k/x/2+alfaratio_t(1)/x)*abs(psi0/x)**2)*x**2, x, lenx, energy)

!!**********************************************************************************************
!***********************************************************************************************


!Prepare Evolution

psi_k=psi0

n=0

!***********************************************************************************************
!*********************************** Perform time iteration ************************************
!***********************************************************************************************
do it=1, lent-1

	if(timearray(it)>tf) then
		goto 444
	endif
	
	call cpu_time(start)

	write(3,*)
	write(2,*)


	pot_r = alfaratio_t(it)/x


	dt = timearray(it+1)-timearray(it)

	timediag = [(0.0d0+img/dt, i=1,lenx)]

	!####################################################################################################
	!################################## Predictor Step Scheme ###########################################
	!####################################################################################################

	call predictor_Article()

	!####################################################################################################
	!####################################################################################################



	!####################################################################################################
	!##################################  Modulation Scheme  #############################################
	!####################################################################################################

		! Modulate the wave function #the self-consistent method does not require this step
		!call Modulation()
	
		! Refine the solution through a self-consistency cycle
		call self_consistent()


		If (omgout == 99999999) then
			Go to 444
		Endif

	!####################################################################################################
	!####################################################################################################

		write(6666,*) timearray(it)
		write(8,*) energy
		write(9,*) norm

		! Write the output values as needed. 

		frac=0.01
		call writer(0)

		! full == 0 --> writes only the instants where alpha ratio are frac=0.01 apart (can be changed)
		!	      These instant are WRITTEN in (66)	(as "tout.txt"), not to be mistaken by 
		!	      (6666) (as "instants.txt") which are ALL the EVALUATED instants.	

		! full != 0 --> writes every instant (frac will not matter). And "instants.txt" should be equiv. to "tout.txt".

		! Note that the more you write the longer it takes to run. Furthermore, because these parameters are not defined 
		! as input (which easily can) any changes require recompiling.


		call cpu_time(finish)

		timetest(it)=(finish-start)

		! Print expected runtime of the simulation.

		call time_verb(0)
		! if input == 0 prints in terminal, else prints in external file (232) (as "TimeLeft.txt").

	enddo

! In case you wish to retrieve the last computed simulation

do i=1,lenx
	write(4321,*) psi_k(i)
enddo

444 END PROGRAM CrankNicolson


