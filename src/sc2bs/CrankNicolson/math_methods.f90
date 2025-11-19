subroutine solve_tridiag(a,b,c,d,x,n)
      implicit none
!	 a - sub-diagonal 
!	 b - the main diagonal
!	 c - sup-diagonal 
!	 d - right part
!	 n - number of equations
!	 x - solution

        integer,intent(in) :: n
	real*8, dimension(n), intent(in) :: a,c
        complex*16,dimension(n),intent(in) :: b,d
        complex*16,dimension(n),intent(out) :: x
        complex*16,dimension(n) :: cp,dp
        complex*16 :: m
        integer:: i


! initialize c-prime and d-prime
        cp(1) = c(1)/b(1)
        dp(1) = d(1)/b(1)
	
! solve for vectors c-prime and d-prime
         do i = 2,n
           m = b(i)-cp(i-1)*a(i)
           cp(i) = c(i)/m
           dp(i) = (d(i)-dp(i-1)*a(i))/m
         end do

! initialize x
         x(n) = dp(n)
	
! solve for x from the vectors c-prime and d-prime
        do i = n-1, 1, -1
          x(i) = dp(i)-cp(i)*x(i+1)
        end do

end subroutine solve_tridiag

subroutine solve_tridiagr(a, b, c, d, x, v, n)
      implicit none
!	 a - sub-diagonal 
!	 b - main diagonal
!	 c - sup-diagonal 
!	 d - right part
!	 n - number of equations
!	 v - some value for the initial conditions
!	 x - solution

        integer,intent(in) :: n
     	real*8, dimension(n), intent(in) :: a,c,b,d
       
        real*8,dimension(n),intent(out) :: x
        real*8,dimension(n) :: cp,dp
        real*8 :: m,v
        integer:: i


! initialize c-prime and d-prime
        cp(1) = c(1)/b(1)
        dp(1) = d(1)/b(1)
	
! solve for vectors c-prime and d-prime
         do i = 2,n
           m = b(i)-cp(i-1)*a(i)
           cp(i) = c(i)/m
           dp(i) = (d(i)-dp(i-1)*a(i))/m
         end do

! initialize x
! Boundary conditions for the potential
         x(n) = 1 !dp(n)
	 x(n-1) = 1 - v
! solve for x from the vectors c-prime and d-prime
        do i = n-1, 1, -1
          x(i) = dp(i)-cp(i)*x(i+1)
        end do

end subroutine solve_tridiagr

!*************************************************************************/
subroutine integrate(f, x, n, integral)
!*************************************************************************/
	implicit none
	integer, intent(in) :: n
	real*8, dimension(n), intent(in):: f, x
	real*8, intent(out):: integral
	integer :: i,m,j

	real*8 :: h0,h1,hph,hdh,hmh, result
	if(mod(n,2)/=0) then
	write(*,*) "Number of grid points is not even", n
	endif
	integral=0
	do i = 2, n-2, 2
            h0  = x(i)-x(i-1)
            h1  = x(i+1)-x(i)
            hph = h1 + h0      
            hdh = h1 / h0
            hmh = h1 * h0

             integral = integral + (hph / 6.0d0) * ((2.00d0 - hdh)*f(i-1)  + (hph**2 / hmh)*f(i)+ (2.0d0 - 1.0d0/hdh)*f(i+1))
	enddo
	!do i=2, N-1
	!integral=integral+ (f(i+1)+f(i))/2*(x(i+1)-x(i))
	!if(mod(i,2)==0) then
	!	integral=integral+ 4.0d0/3.0d0*f(i)
	!else 
	!integral=integral+ 2.0d0/3.0d0*f(i)
	!endif
	
!	enddo
!	integral=(integral+ 1.0d0/3.0d0*(f(1)+f(N)))*(x(3)-x(2))

end subroutine integrate



!*************************************************************************/
subroutine integrateC(f, x, n, integral)
!*************************************************************************/
	implicit none
	integer, intent(in) :: n
	Complex*16, dimension(n), intent(in):: f
	real*8, dimension(n), intent(in):: x
	real*8, intent(out):: integral
	integer :: i,m,j
	real*8 ::   h0,h1,hph,hdh,hmh

	if(mod(n,2)/=0) then
	write(*,*) "Number of grid points is not even", n
	endif


	integral=0
	!do i = 2, n, 2
         !   h0  = x(i)-x(i-1)
          !  h1  = x(i+1)-x(i)
           ! hph = h1 + h0      
            !hdh = h1 / h0
            !hmh = h1 * h0

            !integral = integral + (hph / 6.0d0) * ((2.0d0 - hdh)*f(i-1)       &
             !             + (hph**2 / hmh)*f(i)                            &
             !             + (2.00d0 - 1.00d0/hdh)*f(i+1))
        !end do

	

	do i=2, N-1
	!integral=integral+ (f(i+1)+f(i))/2*(x(i+1)-x(i))
	if(mod(i,2)==0) then
		integral=integral+ 4.0d0/3.0d0*f(i)
	else 
	integral=integral+ 2.0d0/3.0d0*f(i)
	endif
	
	enddo
	integral=(integral+ 1.0d0/3.0d0*(f(1)+f(N)))*(x(3)-x(2))
end subroutine integrateC



!*************************************************************************/
  SUBROUTINE grad(gg, wf, x, llenx)
!**************************************************************************/

Implicit none
	real*8, dimension(llenx), intent(in) :: wf, x
	real*8 , dimension(llenx),intent(out):: gg
        Integer, intent(in):: llenx 
	integer::i
	real*8 :: a,b,c

	do i=2, llenx-1
		a = ( 2*x(i) - x( i ) - x(i+1) ) / ( (x(i-1) - x( i )) * (x(i-1) - x(i+1)))
		b = ( 2*x(i) - x(i-1) - x(i+1) ) / ( (x( i ) - x(i-1)) * (x(i)   - x(i+1)))
		c = ( 2*x(i) - x(i-1) - x( i ) ) / ( (x(i+1) - x(i-1)) * (x(i+1) - x( i )))


		gg(i) =  a*wf(i-1) + b*wf(i) +c*wf(i+1) 
	enddo
	gg(1)=gg(2)
        gg(llenx)=gg(llenx-1)


End Subroutine grad

!*************************************************************************/
  SUBROUTINE gradc(gg, wf, x, llenx)
!**************************************************************************/

Implicit none
	Complex*16, dimension(llenx), intent(in) :: wf 
	real*8, dimension(llenx), intent(in) :: x
	complex*16 , dimension(llenx),intent(out):: gg
        Integer, intent(in):: llenx 
	integer::i
	real*8:: a,b,c

	do i=3, llenx-2
		a = ( 2*x(i) - x( i ) - x(i+1) ) / ( (x(i-1) - x( i )) * (x(i-1) - x(i+1)))
		b = ( 2*x(i) - x(i-1) - x(i+1) ) / ( (x( i ) - x(i-1)) * (x(i)   - x(i+1)))
		c = ( 2*x(i) - x(i-1) - x( i ) ) / ( (x(i+1) - x(i-1)) * (x(i+1) - x( i )))
		!gg(i)=  a*wf(i-1) + b*wf(i) +c*wf(i+1) 

		!gg(i) = wf(i+1)-wf(i)/(x(i+1)-x(i))
		gg(i) = (-wf(i+2)+8*wf(i+1)-wf(i-1)*8+wf(i-2))/(12*( x(i+1)-x(i) ))
	enddo
	
	gg(1)=wf(1)/x(1)
        gg(llenx)=gg(llenx-1)

End Subroutine gradc

!*************************************************************************/
  SUBROUTINE derivative(dfdx, xp, xb, xi, xf, fb, fi, ff)
!**************************************************************************/

!derivative of f at xp

Implicit none
	
	real*8,  intent(in) :: xp, xb, xi, xf, fb,fi,ff
	real*8, intent(out):: dfdx
     
	real*8:: a,b,c
	

	a = ( 2*xp - xi - xf ) / ( (xb - xi) * (xb - xf) )
	b = ( 2*xp - xb - xf ) / ( (xi - xb) * (xi - xf) )
	c = ( 2*xp - xb - xi ) / ( (xf - xb) * (xf - xi) )

	dfdx =  a*fb + b*fi +c*ff

End Subroutine derivative 

!*************************************************************************/
  SUBROUTINE tridiag_mult(s, d, t, llenx, x, vector)
!**************************************************************************/

Implicit none
!	 A.x = vector
!	 s - sub-diagonal 
!	 d - the main diagonal
!	 t - top-diagonal 
	integer, intent(in) :: llenx
	real*8, dimension(llenx), intent(in) :: s, t
	Complex*16, dimension(llenx), intent(in) :: x, d !, s,t
	complex*16 , dimension(llenx),intent(out):: vector
       	integer:: i
	
	do i=2, llenx-1
		vector(i)= s(i)*x(i-1) + d(i)*x(i) + t(i)*x(i+1) 
	enddo
	!boundary conditions 
	vector(1)= d(1)*x(1)+t(1)*x(2)
	vector(llenx) = s(llenx)*x(llenx-1)+d(llenx)*x(llenx) 
	
End Subroutine tridiag_mult

!*************************************************************************/
  SUBROUTINE tridiag_mult_real(s, d, t, llenx, x, vector)
!**************************************************************************/

Implicit none
!	 A.x = vector
!	 s - sub-diagonal 
!	 d - the main diagonal
!	 t - sup-diagonal 
	integer, intent(in) :: llenx
	real*8, dimension(llenx), intent(in) :: s, t, d
	Complex*16, dimension(llenx), intent(in) :: x
	complex*16 , dimension(llenx),intent(out):: vector
       	integer:: i
	
	do i=2, llenx-1
		vector(i)= s(i)*x(i-1) + d(i)*x(i) + t(i)*x(i+1) 
	enddo
	!boundary conditions 
	vector(1) = 0!d(1)*x(1)+t(1)*x(2)
	vector(llenx) = 0!s(llenx)*x(llenx-1)+d(llenx)*x(llenx) 
	
End Subroutine tridiag_mult_real

!*************************************************************************/
 SUBROUTINE romberg_integrate_uneven(f, x, n, tol, maxiter, result)
!*************************************************************************/
    implicit none
	integer, intent(in)::n
    real*8,  dimension(n) ,intent(in):: x, f
    real*8, intent(out):: result
    real*8, intent(in):: tol
    integer, intent(in)  :: maxiter
    real*8, allocatable :: R(:,:)
    integer*8 :: i, j, m
    real*8 :: sum, err_est, a

    allocate(R(maxiter, maxiter))

    ! Step 1: Initialize the trapezoidal rule for the entire uneven grid
    call trapezoidal_uneven(x, f, n, a)
    R(1, 1) = a

    do i = 2, maxiter
        
        !m = 2**(maxiter - i + 1) + 1
	m = max(2**(i-1) + 1,n)
	
	call trapezoidal_subgrid(x, f, n, m, a)
        R(i, 1) = a
	
        ! Perform Richardson extrapolation
        do j = 2, i
            R(i, j) = R(i, j-1) + (R(i, j-1) - R(i-1, j-1)) / (4.0d0**(j-1) - 1.0d0)
        end do

        ! Estimate the error
        
	err_est = abs((R(i, i) - R(i-1, i-1))/R(i,i))
        if (err_est <=tol) then
            result = R(i, i)
	!write(*,*) result
            deallocate(R)
            return
        end if
    end do

    ! If tolerance is not met, return the last value
	write(*,*) "NOT achieved"
    result = R(maxiter, maxiter)
    deallocate(R)
end subroutine romberg_integrate_uneven


!*************************************************************************/
 SUBROUTINE trapezoidal_uneven(x, f, n, tt) 
!*************************************************************************/
    implicit none
    integer, intent(in) :: n
    real*8, dimension(n) , intent(in) :: x, f
    real*8, intent(out) :: tt
    integer :: i

    tt = 0.0d0
    do i = 1, n-1
        tt = tt + 0.5d0 * (x(i+1) - x(i)) * (f(i) + f(i+1))
    end do
end subroutine trapezoidal_uneven

! Trapezoidal rule on a subgrid
!*************************************************************************/
 SUBROUTINE trapezoidal_subgrid(x, f, n, m, tt) 
!*************************************************************************/
    implicit none
    integer, intent(in) :: n 
	integer*8, intent(in)::m
    real*8, dimension(n), intent(in) :: x, f
    real*8, intent(out) :: tt
    integer :: i, step

    !h = (x(n) - x(1)) / (m - 1)
    tt = 0.0d0
	tt = 0.0d0
    do i = 1, m-1
        tt = tt + 0.5d0 * (x(i+1) - x(i)) * (f(i) + f(i+1))
    end do

   ! do i = 1, m-1
   !     tt = tt + 0.5d0 * (x(1 + (i-1)*step + step) - x(1 + (i-1)*step)) *(f(1 + !(i-1)*step) + f(1 + (i-1)*step + step))
 !   end do
end subroutine trapezoidal_subgrid

!*************************************************************************/
  SUBROUTINE interpolate(f0,x0,f1,x1,f2,x2,xx, y) 
!*************************************************************************/
real*8, intent(in) :: f0,x0,f1,x1,f2,x2,xx
real*8, intent(out) :: y

y = f0*((xx - x1)*(xx - x2))/((x0 - x1)*(x0 - x2)) + f1*((xx - x0)*(xx - x2))/((x1 - x0)*(x1 - x2)) + f2*((xx - x0)*(xx - x1))/((x2 - x0)*(x2 - x1))

end subroutine interpolate

!*************************************************************************/
  SUBROUTINE integral_RK4( f, x, n, integral) 
!*************************************************************************/

real*8 :: f0,x0,f1,x1,f2,x2,y,k1,k2,k3,k4,h
integer :: n
real*8, dimension(n), intent(in):: x,f
real*8, intent(out):: integral
!since f is supposed to be an array for rk4 we need to interpolate

y=0
do i=1, n-2
	x0 = x(i);   x1 = x(i+1);   x2 = x(i+2)
	f0 = f(i);   f1 = f(i+1);   f2 = f(i+2)
	h= x(i+1)-x(i)
	call interpolate(f0, x0, f1, x1, f2, x2, x0, k1 )
	call interpolate(f0, x0, f1, x1, f2, x2, x0 + h/2, k2)
	call interpolate(f0, x0, f1, x1, f2, x2, x0 + h/2, k3 )
	call interpolate(f0, x0, f1, x1, f2, x2, x0 + h, k4)
	y = y + h/6*(k1+ 2*k2 + 2*k3 + k4)

enddo

integral = y


end subroutine integral_rk4

!*************************************************************************/
  SUBROUTINE integral_RK5( f, x, n, integral) 
!*************************************************************************/
real*8 :: f0,x0,f1,x1,f2,x2,y,k1,k2,k3,k4,h,k5, k6
integer :: n
real*8, dimension(n), intent(in):: x,f
real*8, intent(out):: integral
!since f is supposed to be an array for rk4 we need to interpolate

y=0
do i=1, n-2
	x0 = x(i);   x1 = x(i+1);   x2 = x(i+2)
	f0 = f(i);   f1 = f(i+1);   f2 = f(i+2)
	h= x(i+1)-x(i)
	call interpolate(f0, x0, f1, x1, f2, x2, x0, k1 )
	call interpolate(f0, x0, f1, x1, f2, x2+1.0d0/5.0d0*h, x0, k2 )
	call interpolate(f0, x0, f1, x1, f2, x2, x0 + 3.0d0/10.0d0*h, k3)
	call interpolate(f0, x0, f1, x1, f2, x2, x0 + 4.0d0/5.0d0*h, k4)
	call interpolate(f0, x0, f1, x1, f2, x2, x0 + 8.0d0/9.0d0*h, k5)
	call interpolate(f0, x0, f1, x1, f2, x2, x0 + h, k6)
	c1=35.0d0/384.0d0; c2=0.0; c3=500.0d0/1113.0d0;
	c4=125.0d0/192.0d0; c5=-2187.0d0/6784.0d0; c6=11.0d0/84.0d0
	y = y + h * ( c1*k1 + c2*k2 + c3*k3 + c4*k4 + c5*k5 + c6*k6 )

enddo

integral = y


end subroutine integral_rk5

!*************************************************************************/
  SUBROUTINE integral_RK6( f, x, n, integral) 
!*************************************************************************/
real*8 :: f0,x0,f1,x1,f2,x2,y,k1,k2,k3,k4,h,k5,k6
integer :: n
real*8, dimension(n), intent(in):: x,f
real*8, intent(out):: integral
!since f is supposed to be an array for rk4 we need to interpolate

y=0
do i=1, n-2
	x0 = x(i);   x1 = x(i+1);   x2 = x(i+2)
	f0 = f(i);   f1 = f(i+1);   f2 = f(i+2)
	h= x(i+1)-x(i)
	call interpolate(f0, x0, f1, x1, f2, x2, x0, k1 )
	call interpolate(f0, x0, f1, x1, f2, x2, x0 + h/3, k2)
	call interpolate(f0, x0, f1, x1, f2, x2, x0 + h/3, k3 )
	call interpolate(f0, x0, f1, x1, f2, x2, x0 + h/2, k4)
	call interpolate(f0, x0, f1, x1, f2, x2, x0 + h, k5)
	call interpolate(f0, x0, f1, x1, f2, x2, x0 + h, k6)
	y = y + h*(k1/12.0d0+ 0*k2 + 3.0d0/8.0d0*k3 + 3.0d0/8.0d0*k4+0*k5 +k6/12.0d0)

enddo

integral = y


end subroutine integral_rk6

