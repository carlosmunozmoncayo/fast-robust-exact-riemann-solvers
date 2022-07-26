! Exact RS SWEs using Van Leer's method with TS initial guess
! =========================================================
subroutine pstar(n_data, rp_data,conv_criteria, tol, gamma, exec_time)
! =========================================================

    implicit none

    integer, intent(in) :: n_data, conv_criteria
    double precision, dimension(n_data,6), intent(in) :: rp_data
    double precision, intent(in) :: tol, gamma
    double precision, intent(out) :: exec_time

    double precision :: pStarNonVacuum, pStarNonVacuum_stagnation
    double precision :: a_l, a_r
    integer :: iterations, i

    double precision :: start_time, finish_time, dummy_p

    iterations=0

    call cpu_time(start_time)
    do i=1,n_data
        a_l = dsqrt(gamma*rp_data(i,1)/rp_data(i,3))
        a_r = dsqrt(gamma*rp_data(i,4)/rp_data(i,6))
        dummy_p=pStarNonVacuum(gamma,rp_data(i,1),rp_data(i,2), &
                    rp_data(i,3),rp_data(i,4),rp_data(i,5),rp_data(i,6),a_l,a_r,tol,iterations)
    end do 
    call cpu_time(finish_time)
    exec_time=finish_time-start_time

    print *,"Iterations Van Leer TS: ", iterations/(0.d0+n_data)

end subroutine pstar

double precision function TwoRarefactionInitialGuess(gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol) 
  !This the exact solution based on having two rarefactions
  implicit double precision (a-h,o-z)
  gamma1=gamma-1.d0
  TwoRarefactionInitialGuess = ((aL+aR-0.5d0*gamma1*(uR-uL))/(aL/pL**(gamma1/(2*gamma))  &
                               +aR/pR**(gamma1/(2*gamma))))**(2*gamma/gamma1)

END function TwoRarefactionInitialGuess


double precision function pStarNonVacuum(gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol,iterations) 
    !This function consider the case where pressure is positive in the initial states and its positivity is preserved  in the middle and gives an estimate on hStar 
    implicit none

    double precision, intent(in) :: gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol
    integer :: iterations

    double precision :: CL,CR,BL,BR
    double precision ::pMin, phiMin, TwoRarefactionInitialGuess, phi, phip, extra, ustar

    !To count iterations required for a singe Riemann Problem
    integer:: iterations_single_RP
    !To decide whether we iterate or not
    integer:: iterate

    !Our estimate for pstar
    double precision :: p, phiR
    double precision :: ppv, gs

    double precision :: usLp, usRp, ustarL, ustarR

    double precision :: ustarL_prime, ustarR_prime, AB

    !We use CL and CR instead of AL and AR to avoid confusion with sound speeds aL and aR
    CL=2/((gamma+1)*rhoL)
    CR=2/((gamma+1)*rhoR)
    BL=(gamma-1)*pL/(gamma+1)
    BR=(gamma-1)*pR/(gamma+1)

    pMin = min(pL,pR)
    phiMin = phi(gamma,pMin,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)

    iterate = 1 
    iterations_single_RP= 0

    if (0<=phiMin) then 
        ! In this case both waves are rarefactions. We know the solution in this case.
        pStarNonVacuum = TwoRarefactionInitialGuess(gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol)
        iterate = 0
    else
        ! Here we have one rarefaction and one shock or two shocks
        ppv=max(pMin,0.5d0*(pL+pR)-0.125d0*(uR-uL)*(rhoL+rhoR)*(aL+aR))

        p = (gs(gamma, ppv, CL, BL)*pL+gs(gamma, ppv, CR, BR)*pR-(uR-uL))/&
          (gs(gamma,ppv,CL,BL)+gs(gamma,ppv,CR,BR)) 
    END if

    if (iterate==1) then 

        do while(.true.)

            ustarL=uL-(p-pL)/AB(gamma,p,pL,rhoL)
            ustarR=uR+(p-pR)/AB(gamma,p,pR,rhoR)
            usLp=ustarL_prime(gamma, p, pL,rhoL,uL)
            usRp=ustarR_prime(gamma, p, pR,rhoR,uR)

            p=p-(ustarL-ustarR)/(usLp-usRp)

            iterations=iterations+1
            iterations_single_RP=iterations_single_RP+1

            if (abs(1-ustarL/ustarR).LE. tol) then
                exit
             end if

            if (isnan(p)) then 
                print *, "aborted at", iterations_single_RP, "iterations by Nan"
                print *, pL,uL,rhoL,pR,uR,rhoR
                print*, "pstar reached: ", p
                print*, "phiR, phipR:", phiR, phip(gamma,p,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)
                call abort
            END if
            !Check if I have reached the tolerance 


            if (iterations_single_RP>40) then
                print *, "aborted at", iterations_single_RP, "iterations"
                print *, pL,uL,rhoL,pR,uR,rhoR
                print*, "pstar reached: ", p
                print*, "phiR, phipR:", phiR, phip(gamma,p,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)
                call abort
                
            end if

        END do
     pStarNonVacuum = p
    END if
   extra=uStar(gamma,p,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)
END function pStarNonVacuum

double precision function gs(gamma, p, CS, BS)
    implicit none 
    double precision:: gamma, p, CS, BS
    gs=dsqrt(CS/(p+BS))
end function gs

double precision function AB(gamma,pstar,pZ,rhoZ)
    implicit none
    double precision :: gamma,pstar,pZ,rhoZ

    if (pstar .GE. pZ) then
        AB=dsqrt(gamma*rhoZ*pZ)*dsqrt((gamma+1)*pstar/(2*gamma*pZ)+(gamma-1)/(2*gamma))
    else
        AB=((gamma-1)/(2*gamma))*dsqrt(gamma*rhoZ*pZ)*(1-pstar/pZ)/(1-(pstar/pZ)**((gamma-1)/(2*gamma)))
    end if
end function AB

double precision function ustarL_prime(gamma, pstar, pL,rhoL,uL)
    implicit none
    double precision :: a, gamma, pstar, pL,rhoL,uL, AB

   if (pstar .GE. pL) then
        a=AB(gamma,pstar,pL,rhoL)
        ustarL_prime=-(a**2+gamma*rhoL*pL)/(2*a**3)
   else
        ustarL_prime=-(pstar/pL)**((-gamma-1)/(2*gamma))/dsqrt(gamma*rhoL*pL)
   end if
end function

double precision function ustarR_prime(gamma, pstar, pR,rhoR,uR)
    implicit none
    double precision :: b, gamma, pstar, pR,rhoR,uR, AB
    if (pstar .GE. pR) then
        b=AB(gamma,pstar,pR,rhoR)
        ustarR_prime=(b**2+gamma*rhoR*pR)/(2*b**3)
    else
        ustarR_prime=(pstar/pR)**((-gamma-1)/(2*gamma))/dsqrt(gamma*rhoR*pR)
    end if
end function

double precision function uStar(gamma,p,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR) 
    !This is based on TORO(2008) Sec 4.2, Prop 4.1
    !Here we assume that pStar has been found 
    implicit none
    double precision :: gamma,p,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR
    double precision ::fL,fR

    ! Compute fL(p,WL)
    if (p<=pL) then !Rarefaction 
        fL = (2*aL/(gamma-1))*((p/pL)**((gamma-1)/(2*gamma))-1)
    else 
        fL = (p-pL)*dsqrt(CL/(p+BL))
    END if
    ! Compute fR(p,WR) 
    if (p<=pR) then 
        fR = (2*aR/(gamma-1))*((p/pR)**((gamma-1)/(2*gamma))-1)
    else 
        fR = (p-pR)*dsqrt(CR/(p+BR))
    END if

    uStar=(uL+uR+fR-fL)/2.d0
END function uStar

double precision function phi(gamma,p,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)  
    !Pressure function phi(p) = f(p,WL) + f(p,WR) + uR - uL
    !See Toro(2008) Sec 4.2, Prop. 4.1
    implicit none
    double precision :: gamma,p,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR
    double precision :: fL, fR

    ! Compute fL(p,WL)
    if (p<=pL) then !Rarefaction 
        fL = (2*aL/(gamma-1))*((p/pL)**((gamma-1)/(2*gamma))-1)
    else 
        fL = (p-pL)*dsqrt(CL/(p+BL))
    END if
    ! Compute fR(p,WR) 
    if (p<=pR) then 
        fR = (2*aR/(gamma-1))*((p/pR)**((gamma-1)/(2*gamma))-1)
    else 
        fR = (p-pR)*dsqrt(CR/(p+BR))
    END if
    phi = fL + fR + uR - uL
END function phi

double precision function phip(gamma,p,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)
    !Derivative of depth function phi'(p) = fL'(p,WL) + fR'(p,WR)
    implicit none
    double precision :: gamma,p,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR
    double precision :: fpL, fpR

    ! Compute f'(p,WL)
    if (p<=pL) then 
        fpL = (p/pL)**(-(gamma+1)/(2*gamma))/(rhoL*aL) 
    else 
        fpL = dsqrt(CL/(BL+p))*(1-(p-pL)/(2.d0*(BL+p)))
    END if
    ! Compute f'(p,WR) 
    if (p<=pR) then 
        fpR = (p/pR)**(-(gamma+1)/(2*gamma))/(rhoR*aR) 
    else 
        fpR = dsqrt(CR/(BR+p))*(1-(p-pR)/(2.d0*(BR+p)))
    END if
    phip = fpL + fpR  
END function phip