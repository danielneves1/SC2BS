!***************************************************
   MODULE Variables
!***************************************************
Real*8 , Parameter::  mp=2.176E-8, kgtoev=5.609E35,Msun = 1.91691E30
Real*16 , Parameter::GeV = 1.0E9 / kgtoev / mp 
character(len=80) :: output_line  
CHARACTER*150 :: Ficheiro, potout,wfout,xpoints,tout,Energ ,norma, wfin,xout, testsout
REAL*8 , DIMENSION (:), ALLOCATABLE :: V_k, V0, v_k1, v_k2, v_k3,g ,pot_r, grd, dxx
REAL*8 , DIMENSION (:,:), ALLOCATABLE ::   avg_v, D2
REAL*8 , DIMENSION (:), ALLOCATABLE :: x, timearray,F, test,timetest,alfaratio_t,diag,subdiag,topdiag
Complex*16 , DIMENSION (:), ALLOCATABLE :: psi0,psi_k1,  psi_k2, psi_k3, psi_k, work, timediag, wf_grad, Bpsi,wf_grad1
complex*16 , DIMENSION (:,:), ALLOCATABLE :: B, B1, B2

real :: start, finish,factor
INTEGER :: i, j,  lenx, lent, it, info,Nmax, dummy,dummy3, lenx2,n,m
CHARACTER*50 :: spacings, timespacing
Real*8  :: alfaratio, initial_alfa, alfa_c, tau,dx, dt, Mbh, spin, mu, numerator, denominator,omg, perc,dum,tf,xi,xf,omgout,to,res,res1, loss, frac
complex*16, parameter :: img = (0.0, 1.0)
Real*8  ::xn, xnew,slope,dummy2, norm,k,energy,energk, energy_2
complex*16 :: dem, nume

Complex*16, dimension(:), allocatable :: psi_test, psi_testold
real*8, dimension(:), allocatable :: V_test, V_testold
END MODULE Variables