
!**************************************************************************************
Subroutine Modulation()
!**************************************************************************************
Use variables
Implicit none


Call integrate((v_k2/x - v_k1/x) * (abs(psi_k2)**2.0 - abs(psi_k)**2.0 )/dt , x, lenx, numerator )

call gradc(wf_grad,    psi_k2, x, lenx)
call gradc(wf_grad1, g*psi_k2, x, lenx)

Call integrate(4* dt**2.0 * Real((-img)*conjg(wf_grad1) * wf_grad ), x, lenx, denominator )

omg = -numerator/denominator


If (omgout== 99999999) then
write(*,*) "Newton did not converge"
return
Endif

omgout=omg

psi_k = psi_k2 * exp( img* omgout* dt**3 * g)
v_k = v_k2



end Subroutine Modulation




!**************************************************************************************
Subroutine self_consistent()
!**************************************************************************************
Use variables
Implicit none


integer:: maxiter, iters
real*8:: energy_test, tolerance, norm_test, energy_loss, potdiff, energy_test1, idk

!allocate(psi_test(lenx), psi_testold(lenx), V_test(lenx))


psi_test = psi_k2
V_test = V_k2
maxiter=500


do iters=1, maxiter 

	psi_testold = psi_test
	V_testold = V_test
	

	call wf(psi_test, (V_test + V_k + alfaratio_t(it+1) + alfaratio_t(it))/x/2)
	
	call CloudPotential(V_test, psi_test)
	
	call integrate(abs(psi_test)**2, x, lenx, norm_test)
	psi_test = psi_test/ sqrt(norm_test)

	call gradc(wf_grad, psi_test/x, x, lenx)
	
	!Convergence Condition Evaluation
	call integrate(abs(psi_test-psi_testold), x, lenx, tolerance)
	
	
	
!Make it so that at least one refinement is done
!Convergence Condition

if ( iters>1.and.tolerance<1e-12) then

	psi_k = psi_test
	V_k = V_test
	call integrate((abs(wf_grad)**2 - (V_test/2 + alfaratio_t(it+1)) /x * abs(psi_test/x)**2.0)*x**2, x, lenx, energy)
	!energy = energy_test
	norm = norm_test
	return
endif



if(iters==maxiter) then
omgout = 99999999
write(*,*) "Self-Consistent method did not converge :", tolerance
endif
enddo




end Subroutine self_consistent

