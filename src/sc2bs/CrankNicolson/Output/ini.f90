!*************************************************************************/
  SUBROUTINE Initialize
!**************************************************************************/


!For this particular code we need an initial wave-function alongside the spatial grid (x-grid) where it is defined.
!Inside the x-grid the first line should only have an integer corresponding to the total length of the x-grid( and thus the wf)


Use Variables


IMPLICIT NONE

!************************************************************************************************
!***************************** Read Input Paramenters *******************************************
!************************************************************************************************

!Most x-grid of these parameters are legacy parameters and can be used if knowledge of the initial wave function
! is known analytically. In this way, we can define the initial wave-function here for the input x-grid from defined 
! by starting from xi to xf with lenx points, with the option of having factor * lenx of these points written in log-scale.
! That is, if defined.

!*******************************************************************************************************
!The real important parameters for this file are the time parameters and alpha ratio and tau parameters
!*******************************************************************************************************

READ (1,*) spacings
READ (1,*) !lenx no need since an x file must be provided
READ (1,*) !xi !NOT USED
READ (1,*) !xf !NOT USED
READ (1,*) factor

READ(1,*)
READ(1,*)  timespacing
Read(1,*)  Nmax
Read(1,*)  lent
Read(1,*)  to
Read(1,*)  tf	
READ(1,*)  k !minimum variation of the Black Hole-Cloud mass ratio (see SUBROUTINE tempos in instants.f90 file )

READ(1,*)  alfaratio
Read(1,*)
Read(1,*) Mbh !NOT USED
Read(1,*) spin !NOT USED
Read(1,*) tau
Read(1,*)
Read(1,*) mu !NOT USED



read(4,*) lenx

!************************************************************************************************
!************************************************************************************************




!************************************************************************************************
!***************************** Define Spatial Grid Quantities ***********************************
!************************************************************************************************

if( mod(lenx,2) /=0 ) then
	lenx=lenx+1
end if

allocate(x(lenx))

do i=1, lenx-1

	Read(4,*) x(i)  
enddo

x(lenx) = x(lenx-1) + x(lenx-1)-x(lenx-2)


do i=1, lenx
	write(22,*) x(i)  
enddo

allocate(dxx(lenx))

do i=1, lenx-1
	dxx(i)=x(i+1)-x(i)
enddo

dxx(lenx)=dxx(lenx-1)


allocate(V_k1(lenx))
allocate(V_k2(lenx))
allocate(V_k(lenx), psi_k(lenx))
allocate(psi_k1(lenx),psi_k2(lenx), wf_grad(lenx))
allocate(pot_r(lenx))
allocate(grd(lenx))
allocate(subdiag(lenx),topdiag(lenx),diag(lenx))
allocate(psi0(lenx))
allocate(Bpsi(lenx), wf_grad1(lenx))


!If memory is a concern, then it should noted that not all allocatables 
! are required. Here they are only defined for better/clearer understanding of the numerical procedure.



!***********************************************************************************
!******************** Define the Finite Differences Method arrays ******************
!***********************************************************************************
do i=2, lenx-1

	subdiag(i) = 2/((x(i-1) - x(i  )) * (x(i-1) - x(i+1)))
	diag(i) =    2/((x( i ) - x(i-1)) * (x(i  ) - x(i+1)))
	topdiag(i) = 2/((x(i+1) - x(i-1)) * (x(i+1) - x(i  )))
enddo

diag(1)=1.0d0
diag(lenx)=1.0d0
subdiag(lenx)=0.0d0
subdiag(1)=0.0d0
topdiag(1)=0.0d0
topdiag(lenx)=0.0d0

!***********************************************************************************
!This code is set for even grid-spacing but somewhat prepared for uneven-grid spacing
!***********************************************************************************


!Reading initial wave function

do i=1, lenx-1
	read(70 ,*) psi0(i)
enddo
psi0(lenx)=0.0*img


call integrate(abs(psi0)**2, x, lenx, norm)


write(232,*) "Initial norm = ", norm

g = 1/(1+x)
!g= Cos(x)


call grad(grd, g, x, lenx)


!************************************************************************************************
!************************************************************************************************


!************************************************************************************************
!***************************** Define Time Dependent Quantities *********************************
!************************************************************************************************


!This computes the required instant such that the BlackHole-Cloud mass ratio evolves more smoothly.
!This means that the timespacing is uneven which should always be the case for sharp variations of
! the potential. This is achieved by calling the SUBROUTINE "tempos" which should be modified for 
! different evolution laws.

if(timespacing=="UnEven") then
	allocate(timearray(Nmax+lent))
	timearray(1)=to

	do i=1,Nmax
		call tempos(dum, timearray(i), k)
	
		if( dum/tf>=1) then
			j=i

			goto 101
		endif
	
		timearray(i+1)= to+dum
		if(i==Nmax) then
			write(*,*) "Not enough time steps allowed for k =", k
			write(*,*) "Either consider, in input file, changing Max Num. of time steps before evap =", Nmax, "or change k."

			return
		endif
	enddo

	if(i==Nmax) then
		write(232,*) "Not enough time interval rangevalues"
		return
	endif

else 
	j=0
	allocate(timearray(lent))

	write(232,*)

endif



!After the external potential (Black Hole) evaporates then the time spacing can be even
! and adjusted to capture the wave function evolution. Here the input time length alongside
! the final instant (also input) defines the step.

101 dt = (tf-timearray(j))/lent

do i=1,lent
timearray(j+i)=timearray(j) + i*dt
enddo

write(232,"(A,Es12.4, A)") " Time interval dt  = ",dt, "(after Evaporation)"

lent=j+lent

allocate(alfaratio_t(lent))
allocate(timetest(lent-1))


!************************************************************************************************
!************************************************************************************************
	!Obtain Black Hole - Cloud mass ratio evolution for each evaluated instant 	
!************************************************************************************************
!************************************************************************************************

	do i=1, lent
		if((timearray(i)-to)<tau) then
		alfaratio_t(i) = alfaratio * (1 - (timearray(i)-to)/tau)**(1.0d0/3.0d0)

		else
			alfaratio_t(i)=0.0d0
		endif
	enddo

!************************************************************************************************
!************************************************************************************************


!************************************************************************************************
!************************************************************************************************



write(232,*)

write(232,*) "SpaceGrid Size:", lenx
write(232,*) "TimeGrid Size:", lent


timetest=0.0

End Subroutine Initialize





