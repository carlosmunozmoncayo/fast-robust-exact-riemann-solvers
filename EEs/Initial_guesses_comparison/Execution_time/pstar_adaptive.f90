! Exact RS SWEs using Newton's method with Convex Combination initial guess
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

    if (conv_criteria==1) then
        !We use a convergence criteria based on the residual
        call cpu_time(start_time)
        do i=1,n_data
            a_l = dsqrt(gamma*rp_data(i,1)/rp_data(i,3))
            a_r = dsqrt(gamma*rp_data(i,4)/rp_data(i,6))
            dummy_p=pStarNonVacuum(gamma,rp_data(i,1),rp_data(i,2), &
               rp_data(i,3),rp_data(i,4),rp_data(i,5),rp_data(i,6),a_l,a_r,tol,iterations)
            !print *, dummy_p
        END do 
        call cpu_time(finish_time)
        exec_time=finish_time-start_time
        print *,"Iterations Newton CC: ", iterations/(n_data +0.d0)

    else if (conv_criteria==0 ) then
        !We use a convergence criteria based on stagnation
        call cpu_time(start_time)
        do i=1,n_data
            a_l = dsqrt(gamma*rp_data(i,1)/rp_data(i,3))
            a_r = dsqrt(gamma*rp_data(i,4)/rp_data(i,6))
            dummy_p=pStarNonVacuum_stagnation(gamma,rp_data(i,1),rp_data(i,2), &
               rp_data(i,3),rp_data(i,4),rp_data(i,5),rp_data(i,6),a_l,a_r,tol,iterations)
            !print *, dummy_p
        END do 
        call cpu_time(finish_time)
        exec_time=finish_time-start_time
        print *,"Iterations Newton CC stagnation criteria: ", iterations/(n_data +0.d0)

    else
        print *, "Convergence criteria not properly specified"
        call abort

    END if

end subroutine pstar

double precision function TwoRarefactionInitialGuess(gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol) 
    !This the exact solution based on having two rarefactions
    implicit none
    double precision, intent (in) :: gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol
    double precision :: gamma1
    gamma1=gamma-1.d0
    TwoRarefactionInitialGuess = ((aL+aR-0.5d0*gamma1*(uR-uL))/(aL/pL**(gamma1/(2*gamma)) &
        +aR/pR**(gamma1/(2*gamma))))**(2*gamma/gamma1)

END function TwoRarefactionInitialGuess


double precision function pStarNonVacuum(gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol,iterations) 

    implicit none

    double precision, intent(in) :: gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol
    integer :: iterations

    double precision :: CL,CR,BL,BR
    double precision ::pMin, phiMin, TwoRarefactionInitialGuess, phi, phip, extra, ustar
    double precision ::pMax, phiMax

    !To count iterations required for a singe Riemann Problem
    integer:: iterations_single_RP
    !To decide whether we iterate or not
    integer:: iterate

    !Required for initial guess
    double precision :: ppv, gs

    !Our estimate for pstar
    double precision :: p, phiR


    !We use CL and CR instead of AL and AR to avoid confusion with sound speeds aL and aR
    CL=2/((gamma+1)*rhoL)
    CR=2/((gamma+1)*rhoR)
    BL=(gamma-1)*pL/(gamma+1)
    BR=(gamma-1)*pR/(gamma+1)

    ! We estimate pstar from below, this is useful even if we need to iterate
    pMin = min(pL,pR)
    phiMin = phi(gamma,pMin,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)

    iterate = 1 
    iterations_single_RP= 0

    if (0<=phiMin) then 
        ! In this case both waves are rarefactions. The solution is explicit
        pStarNonVacuum = TwoRarefactionInitialGuess(gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol)
        iterate = 0

    else 
        pMax = max(pL,pR)
        phiMax = phi(gamma,pMax,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)

        if (phiMax < 0) then 
            ! In this case pMax < pStar so both waves are shocks.
            p = TwoRarefactionInitialGuess(gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol) 
            ! this estimate is guaranteed to be from the right 
            phiR=phi(gamma,p,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)
            p =(phiR*pMax-phiMax*p)/(phiR-phiMax)
        else
            ! Here we have one rarefaction and one shock 
            ! A left estimate is pMin and a right estimate is min(pMax,p_TR)
            p = min(pMax,TwoRarefactionInitialGuess(gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol))
            phiR=phi(gamma,p,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)
            pStarNonVacuum =(phiR*pMin-phiMin*p)/(phiR-phiMin) 
            iterate = 0
        END if
    END if

    if (iterate==1) then 
        p = max(pMax,p-phi(gamma,p,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)/&
            phip(gamma,p,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR))     
        iterations=iterations+1

        !Start iterative process 
        do while(.true.)
            
            phiR = phi(gamma,p,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)
            
            !Check if we have reached the tolerance 
            if (abs(phiR)<tol) then 
                exit
            END if

            !Compute new estimates of pStar
            p = p-phiR/&
                phip(gamma,p,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)

            iterations=iterations+1
            iterations_single_RP=iterations_single_RP+1

            if (isnan(p).or.isnan(phiR)) then 
                print *, "aborted at", iterations_single_RP, "iterations by Nan"
                print *, pL,uL,rhoL,pR,uR,rhoR
                print*, "pstar reached: ", p
                print*, "phiR, phipR:", phiR, phip(gamma,p,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)
                call abort
            END if
            

            if (iterations_single_RP>40) then
                print *, "aborted at", iterations_single_RP, "iterations"
                print *, pL,uL,rhoL,pR,uR,rhoR
                print*, "pstar reached: ", p
                print*, "phiR, phipR:", phiR, phip(gamma,p,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)
                call abort  
            END if

        END do
        ! return estimation of pStar
        pStarNonVacuum = p
    END if
    extra=uStar(gamma,p,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)
END function pStarNonVacuum

double precision function pStarNonVacuum_stagnation(gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol,iterations) 

    implicit none

    double precision, intent(in) :: gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol
    integer:: iterations

    double precision :: CL,CR,BL,BR
    double precision ::pMin, phiMin, pMax, phiMax, TwoRarefactionInitialGuess
    double precision :: phi, phip, extra, uStar

    !To count iterations required for a singe Riemann Problem
    integer:: iterations_single_RP
    !To decide whether we iterate or not
    integer:: iterate

    !Our estimate for pstar
    double precision :: p, pold, phiR



    !We use CL and CR instead of AL and AR to avoid confusion with sound speeds aL and aR
    CL=2/((gamma+1)*rhoL)
    CR=2/((gamma+1)*rhoR)
    BL=(gamma-1)*pL/(gamma+1)
    BR=(gamma-1)*pR/(gamma+1)

    ! We need to start with two states of h. One to the left and one to the right of hStar
    pMin = min(pL,pR)
    pMax = max(pL,pR)
    phiMin = phi(gamma,pMin,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)
    phiMax = phi(gamma,pMax,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)

    iterate = 1 
    iterations_single_RP= 0
    if (0<=phiMin) then 
        ! In this case both waves are rarefactions. We know the solution in this case.
        pStarNonVacuum_stagnation = TwoRarefactionInitialGuess(gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol)
        iterate = 0
    else if (phiMax < 0) then 
        ! In this case pMax < pStar so both waves are shocks.
        p = TwoRarefactionInitialGuess(gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol) ! this estimate is guaranteed to be from the right 
        phiR=phi(gamma,p,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)
        p =(phiR*pMax-phiMax*p)/(phiR-phiMax)
    else
        ! Here we have one rarefaction and one shock 
        ! A left estimate is pMin and a right estimate is min(pMax,p_TR)
        p = min(pMax,TwoRarefactionInitialGuess(gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol))
        phiR=phi(gamma,p,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)
        p =(phiR*pMin-phiMin*p)/(phiR-phiMin) 
    END if

    if (iterate==1) then 
        p = max(pMin,p-phi(gamma,p,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)/&
            phip(gamma,p,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR))     
        iterations=iterations+1
        !Start iterative process 
        do while(.true.)

            !Check if pold and p are close enough
            if (abs(2*(pold-p)/(pold+p)) <tol) then
                exit
            end if
            !Save old estimates of pStar 
            pold = p

            !Compute new estimates of pStar
            phiR = phi(gamma,p,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)
            p = p-phiR/&
                phip(gamma,p,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)
            iterations=iterations+1
            iterations_single_RP=iterations_single_RP+1

            if (isnan(p).or.isnan(phiR)) then 
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
            END if

        END do
        !return estimation of pStar
        pStarNonVacuum_stagnation = p
    END if
    extra=uStar(gamma,p,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR) 
END function pStarNonVacuum_stagnation

double precision function gs(gamma, p, CS, BS)
   implicit double precision (a-h,o-z)
   gs=dsqrt(CS/(p+BS))
end function gs

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

    implicit double precision (a-h,o-z)

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
    implicit double precision (a-h,o-z)
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


