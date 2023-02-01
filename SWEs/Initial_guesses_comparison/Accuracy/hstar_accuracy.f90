!Testing accuracy of initial guesses for SWEs
! =========================================================
subroutine hstar(n_data, rp_data, tol, grav, exec_time)
! =========================================================

    implicit none

    integer, intent(in) :: n_data
    double precision, dimension(n_data,4), intent(in) :: rp_data
    double precision, intent(in) :: tol, grav
    double precision, intent(out)  :: exec_time

    double precision :: hStarWetStates, TwoRarefactionInitialGuess
    double precision :: hStarWetStates_stagnation
    integer :: iterations,i

    double precision :: hStarWetStates_AV,hStarWetStates_QA,hStarWetStates_TR
    double precision :: hStarWetStates_CC,hStarWetStates_TS, hStarWetStates_PV
    double precision :: hStarWetStates_HLLE, hStarAdaptiveNoniterative


    double precision :: start_time, finish_time, dummy_h
    double precision :: app_AV, app_CC, app_QA, app_TR, app_TS, app_PV, app_HLLE, app_ANI
    double precision :: err_AV, err_CC, err_QA, err_TR, err_TS, err_PV, err_HLLE, err_ANI
    double precision :: exact
 
    iterations=0
    err_AV=0.d0
    err_QA=0.d0
    err_CC=0.d0
    err_TR=0.d0
    err_PV=0.d0
    err_TS=0.d0
    err_HLLE=0.d0
    err_ANI=0.d0 !Adaptive noniterative error

    !We use a convergence criteria based on the residual
    call cpu_time(start_time)
    do i=1,n_data
        exact=hStarWetStates(grav, rp_data(i,1), rp_data(i,2), &
                rp_data(i,3), rp_data(i,4), tol,iterations)
        app_AV=hStarWetStates_AV(grav, rp_data(i,1), rp_data(i,2), &
                rp_data(i,3), rp_data(i,4), tol,iterations)
        app_CC=hStarWetStates_CC(grav, rp_data(i,1), rp_data(i,2), &
                rp_data(i,3), rp_data(i,4), tol,iterations)
        app_QA=hStarWetStates_QA(grav, rp_data(i,1), rp_data(i,2), &
                rp_data(i,3), rp_data(i,4), tol,iterations)
        app_TR=hStarWetStates_TR(grav, rp_data(i,1), rp_data(i,2), &
                rp_data(i,3), rp_data(i,4), tol,iterations)
        app_TS=hStarWetStates_TS(grav, rp_data(i,1), rp_data(i,2), &
                rp_data(i,3), rp_data(i,4), tol,iterations)
        app_PV=hStarWetStates_PV(grav, rp_data(i,1), rp_data(i,2), &
                rp_data(i,3), rp_data(i,4), tol,iterations)
        app_HLLE=hStarWetStates_HLLE(grav, rp_data(i,1), rp_data(i,2), &
                rp_data(i,3), rp_data(i,4), tol,iterations)
        
        app_ANI=hStarAdaptiveNoniterative(grav, rp_data(i,1), rp_data(i,2), &
                rp_data(i,3), rp_data(i,4), tol,iterations) !Adaptive noniterative approximation

        err_AV=err_AV+abs(app_AV-exact)/abs(exact)
        err_QA=err_QA+abs(app_QA-exact)/abs(exact)
        err_CC=err_CC+abs(app_CC-exact)/abs(exact)
        err_TR=err_TR+abs(app_TR-exact)/abs(exact)
        err_TS=err_TS+abs(app_TS-exact)/abs(exact)
        err_PV=err_PV+abs(app_PV-exact)/abs(exact)
        err_HLLE=err_HLLE+abs(app_HLLE-exact)/abs(exact)
        
        err_ANI=err_ANI+abs(app_ANI-exact)/abs(exact)
    end do 
    call cpu_time(finish_time)
    exec_time=finish_time-start_time

    print *, "Error AV:", err_AV/(n_data +0.d0)
    print *, "Error QA:", err_QA/(n_data +0.d0)
    print *, "Error CC:", err_CC/(n_data +0.d0)
    print *, "Error TR:", err_TR/(n_data +0.d0)
    print *, "Error TS:", err_TS/(n_data +0.d0)
    print *, "Error PV:", err_PV/(n_data +0.d0)
    print *, "Error HLLE:", err_HLLE/(n_data +0.d0)
    
    print *, "Error ANI:", err_ANI/(n_data +0.d0)


end subroutine hstar

double precision function hStarWetStates(g,hL,hR,uL,uR,tol,iterations)
    !This function consider the case where both states are wet and gives an estimate on hStar 
    implicit none
    integer :: iterate, iterations, iterations_single_RP
    double precision:: g,hL,hR,uL,uR,tol
    double precision:: hMin, hMax, fMin, fMax
    double precision:: h,phi, phip, phiR
    double precision:: TwoRarefactionInitialGuess
    double precision :: yz

    ! We estimate hstar from below
    hMin = min(hL,hR)
    ! hMax = max(hL,hR)
    fMin = phi(g,hMin,hL,hR,uL,uR)
    ! fMax = phi(g,hMax,hL,hR,uL,uR)

    iterate = 1 
    iterations_single_RP=0
    if (0<=fMin) then 
        ! In this case both waves are rarefactions. We know the solution in this case.
        hStarWetStates = TwoRarefactionInitialGuess(g,hL,hR,uL,uR)
        iterate = 0
    else
        h=0.5d0*(hL+hR)
        h=(hL*yz(g,h,hL)+hR*yz(g,h,hR)-uR+uL)/(yz(g,h,hL)+yz(g,h,hR))
    END if


    if (iterate==1) then 
        !We approach the solution through Newton's method
        !h_{i+1}=h_i-phi(h_i)/phi'(h_i)
        h=max(hMin, h-phi(g,h,hL,hR,uL,uR)/phip(g,h,hL,hR))
        iterations_single_RP=iterations_single_RP+1
        iterations=iterations+1
        phiR=phi(g,h,hL,hR,uL,uR)
        !Start iterative process 
        do while(.true.)
            !Trying to fix 
            if (abs(phiR) <tol) then
                exit
            end if
            h = h-phiR/phip(g,h,hL,hR)   !Updating hstar
            phiR=phi(g,h,hL,hR,uL,uR)              !Updating phi(hstar)
            iterations=iterations+1
            iterations_single_RP=iterations_single_RP+1

            if (isnan(h).or. isnan(phiR)) then
                print *, "Aborted because of Nan, hstar reached",h
                print *, "hL, hR, uL, uR: ", hL, hR, uL, uR
                print *, "Phi, phi prime: ", phiR, phip(g,h,hL,hR)
                call abort
            else if (iterations_single_RP>40) then
                print *, "Exceeded  number of iterations, hstar reached:",h
                print *, "hL, hR, uL, uR: ", hL, hR, uL, uR
                print *, "Phi, phi prime: ", phiR, phip(g,h,hL,hR)
                call abort
            end if
        END do
        ! return estimation of hStar
        hStarWetStates = h
    END if  
END function hStarWetStates

double precision function hStarWetStates_AV(g,hL,hR,uL,uR,tol,iterations)
    !This function consider the case where both states are wet and gives an estimate on hStar 
    implicit none
    integer :: iterate, iterations, iterations_single_RP
    double precision:: g,hL,hR,uL,uR,tol
    double precision:: hMin, hMax, fMin, fMax
    double precision:: h,phi, phip, phiR
    double precision:: TwoRarefactionInitialGuess

    ! We estimate hstar from below
    hMin = min(hL,hR)
    !hMax = max(hL,hR)
    fMin = phi(g,hMin,hL,hR,uL,uR)
    !fMax = phi(g,hMax,hL,hR,uL,uR)

    if (0<=fMin) then 
        ! In this case both waves are rarefactions. We know the solution in this case.
        hStarWetStates_AV = TwoRarefactionInitialGuess(g,hL,hR,uL,uR)
    else
        hStarWetStates_AV=0.5d0*(hL+hR)
    END if 
END function hStarWetStates_AV

double precision function hStarWetStates_QA(g,hL,hR,uL,uR,tol,iterations)
    !This function consider the case where both states are wet and gives an estimate on hStar 
    implicit none
    integer :: iterate, iterations, iterations_single_RP
    double precision:: g,hL,hR,uL,uR,tol
    double precision:: hMin, hMax, fMin, fMax
    double precision:: h,phi, phip, phiR
    double precision:: TwoRarefactionInitialGuess, SharpInitialGuess

    ! We estimate hstar from below
    hMin = min(hL,hR)
    !hMax = max(hL,hR)
    fMin = phi(g,hMin,hL,hR,uL,uR)
    !fMax = phi(g,hMax,hL,hR,uL,uR)

    if (0<=fMin) then 
        ! In this case both waves are rarefactions. We know the solution in this case.
        hStarWetStates_QA = TwoRarefactionInitialGuess(g,hL,hR,uL,uR)
    else
        hStarWetStates_QA=SharpInitialGuess(g,hL,hR,uL,uR)
    END if
END function hStarWetStates_QA

double precision function hStarWetStates_TR(g,hL,hR,uL,uR,tol,iterations)
    !This function consider the case where both states are wet and gives an estimate on hStar 
    implicit none
    integer :: iterate, iterations, iterations_single_RP
    double precision:: g,hL,hR,uL,uR,tol
    double precision:: hMin, fMin
    double precision:: h,phi, phip, phiR
    double precision:: TwoRarefactionInitialGuess

    ! We estimate hstar from below
    hMin = min(hL,hR)
    !hMax = max(hL,hR)
    fMin = phi(g,hMin,hL,hR,uL,uR)
    !fMax = phi(g,hMax,hL,hR,uL,uR)

    hStarWetStates_TR=TwoRarefactionInitialGuess(g,hL,hR,uL,uR)

END function hStarWetStates_TR

double precision function hStarWetStates_HLLE(g,hL,hR,uL,uR,tol,iterations)
    !This function consider the case where both states are wet and gives an estimate on hStar 
    implicit none
    integer :: iterate, iterations, iterations_single_RP
    double precision:: g,hL,hR,uL,uR,tol
    double precision:: hMin, fMin
    double precision:: h,phi, phip, phiR
    double precision:: TwoRarefactionInitialGuess
    double precision:: h_hat, u_hat, c_hat, lambda_1_l, lambda_2_r
    double precision:: s1, s2

    ! We estimate hstar from below
    hMin = min(hL,hR)
    !hMax = max(hL,hR)
    fMin = phi(g,hMin,hL,hR,uL,uR)
    !fMax = phi(g,hMax,hL,hR,uL,uR)

    if (0<=fMin) then 
        ! In this case both waves are rarefactions. We know the solution in this case.
        hStarWetStates_HLLE = TwoRarefactionInitialGuess(g,hL,hR,uL,uR)
    else
        h_hat = (hR + hL)/2.d0
        u_hat = (dsqrt(hR)*uR + dsqrt(hL)*uL) / (dsqrt(hR) + dsqrt(hL))
        c_hat = dsqrt(g*h_hat)

        lambda_1_l = uL - dsqrt(g*hL)
        lambda_2_r = uR + dsqrt(g*hR)
        
        s1 = min(lambda_1_l, u_hat - c_hat)
        s2 = max(lambda_2_r, u_hat + c_hat)
        
        hStarWetStates_HLLE = (hR*uR - hL*uL - s2*hR + s1*hL)/(s1-s2)
    END if

END function hStarWetStates_HLLE

double precision function hStarWetStates_CC(g,hL,hR,uL,uR,tol,iterations)
    !This function consider the case where both states are wet and gives an estimate on hStar 
    implicit none
    integer :: iterate, iterations, iterations_single_RP
    double precision:: g,hL,hR,uL,uR,tol
    double precision:: hMin, hMax, fMin, fMax
    double precision:: h,phi, phip, phiR
    double precision:: TwoRarefactionInitialGuess

    ! We estimate hstar from below
    hMin = min(hL,hR)
    hMax = max(hL,hR)
    fMin = phi(g,hMin,hL,hR,uL,uR)
    fMax = phi(g,hMax,hL,hR,uL,uR)

    if (0<=fMin) then 
        ! In this case both waves are rarefactions. We know the solution in this case.
        hStarWetStates_CC = TwoRarefactionInitialGuess(g,hL,hR,uL,uR)

    else if (fMax < 0) then 
        ! In this case hMax < hStar so both waves are shocks. 
        ! We need to iterate. Here two estimates from the left and right are:
        h=TwoRarefactionInitialGuess(g,hL,hR,uL,uR) !upper bound, lower bound is hmax
        phiR=phi(g,h,hL,hR,uL,uR)
        hStarWetStates_CC = (phiR*hMax-fMax*h)/(phiR-fMax)
    else
        ! Here we have one rarefaction and one shock  
        h = min(hMax,TwoRarefactionInitialGuess(g,hL,hR,uL,uR)) !upper bound, lower bound is hmin
        phiR=phi(g,h,hL,hR,uL,uR)
        hStarWetStates_CC = (phiR*hMin-fMin*h)/(phiR-fMin)
    END if

END function hStarWetStates_CC

double precision function hStarWetStates_PV(g,hL,hR,uL,uR,tol,iterations)
    !This function consider the case where both states are wet and gives an estimate on hStar 
    implicit none
    integer :: iterate, iterations, iterations_single_RP
    double precision:: g,hL,hR,uL,uR,tol
    double precision:: hMin, hMax, fMin, fMax
    double precision:: h,phi, phip, phiR
    double precision:: TwoRarefactionInitialGuess


    ! We estimate hstar from below
    hMin = min(hL,hR)
    ! hMax = max(hL,hR)
    fMin = phi(g,hMin,hL,hR,uL,uR)
    ! fMax = phi(g,hMax,hL,hR,uL,uR)

    if (0<=fMin) then 
        ! In this case both waves are rarefactions. We know the solution in this case.
        hStarWetStates_PV= TwoRarefactionInitialGuess(g,hL,hR,uL,uR)
    else
        hStarWetStates_PV=0.5d0*(hL+hR)+0.25d0*(uL-uR)*(hL+hR)/(dsqrt(g*hL)+dsqrt(g*hR))
    END if


END function hStarWetStates_PV

double precision function hStarWetStates_TS(g,hL,hR,uL,uR,tol,iterations)
    !This function consider the case where both states are wet and gives an estimate on hStar 
    implicit none
    integer :: iterate, iterations, iterations_single_RP
    double precision:: hStarWetStates_PV
    double precision:: g,hL,hR,uL,uR,tol
    double precision:: hMin, hMax, fMin, fMax
    double precision:: h,phi, phip, phiR
    double precision:: TwoRarefactionInitialGuess
    double precision :: yz

    ! We estimate hstar from below
    hMin = min(hL,hR)
    ! hMax = max(hL,hR)
    fMin = phi(g,hMin,hL,hR,uL,uR)
    ! fMax = phi(g,hMax,hL,hR,uL,uR)

    if (0<=fMin) then 
        ! In this case both waves are rarefactions. We know the solution in this case.
        hStarWetStates_TS = TwoRarefactionInitialGuess(g,hL,hR,uL,uR)
    else
        h=hStarWetStates_PV(g,hL,hR,uL,uR,tol,iterations)
        hStarWetStates_TS=(hL*yz(g,h,hL)+hR*yz(g,h,hR)-uR+uL)/(yz(g,h,hL)+yz(g,h,hR))
    END if

END function hStarWetStates_TS

double precision function hStarAdaptiveNoniterative(g,hL,hR,uL,uR,tol,iterations)
    !This function consider the case where both states are wet and gives an estimate on hStar 
    implicit none
    integer :: iterate, iterations, iterations_single_RP
    double precision:: g,hL,hR,uL,uR,tol
    double precision:: hMin, hMax, fMin, fMax
    double precision:: h,phi, phip, phiR
    double precision:: TwoRarefactionInitialGuess, SharpInitialGuess

    ! We estimate hstar from below
    hMin = min(hL,hR)
    !hMax = max(hL,hR)
    fMin = phi(g,hMin,hL,hR,uL,uR)
    !fMax = phi(g,hMax,hL,hR,uL,uR)

    if (0<=fMin) then 
        ! In this case both waves are rarefactions. We know the solution in this case.
        hStarAdaptiveNoniterative = TwoRarefactionInitialGuess(g,hL,hR,uL,uR)
    else
        hMax = max(hL,hR)
        fMax = phi(g,hMax,hL,hR,uL,uR)
        if (fMax<0) then
            hStarAdaptiveNoniterative = SharpInitialGuess(g,hL,hR,uL,uR)
        else
            ! Here we have one rarefaction and one shock  
            h = min(hMax,TwoRarefactionInitialGuess(g,hL,hR,uL,uR)) !upper bound, lower bound is hmin
            phiR = phi(g,h,hL,hR,uL,uR)
            hStarAdaptiveNoniterative = (phiR*hMin-fMin*h)/(phiR-fMin)
        END if
    END if
END function hStarAdaptiveNoniterative

double precision function yz(g, h, hz)
  implicit none
  double precision :: g, h, hz
  yz=dsqrt(0.5d0*g*(h+hz)/(h*hz))
end function yz

double precision function TwoRarefactionInitialGuess(g,hL,hR,uL,uR)
    !This the solution based on having two rarefactions
    ! Note: this can be proven to be an upper bound for hStar
    implicit none
    double precision :: g,hL,hR,uL,uR
    TwoRarefactionInitialGuess = (uL-uR+2.d0*sqrt(g)*(sqrt(hL)+sqrt(hR)))**2.d0/16.d0/g
END function TwoRarefactionInitialGuess

double precision function SharpInitialGuess(g,hL,hR,uL,uR)
    !This is a sharper initial guess for hStar
    !It is designed to work for when one state is near to be dry and the other is not
    implicit none
    double precision :: g, hL, hR, uL, uR
    double precision :: hMin, hMax, x0, fMin, fMax
    double precision :: phi, TwoRarefactionInitialGuess

    !Compute h min/max and the corresponding depth functions
    hMin = min(hL,hR)
    hMax = max(hL,hR)
    x0 = (2.d0*sqrt(2.d0)-1)**2.d0
    fMin = phi(g,x0*hMin,hL,hR,uL,uR)
    fMax = phi(g,x0*hMax,hL,hR,uL,uR)
    !Compute initial guess
    if (0<=fMin) then
     SharpInitialGuess = TwoRarefactionInitialGuess(g,hL,hR,uL,uR)
    else if (0<=fMax) then 
     SharpInitialGuess = (-sqrt(2.d0*hMin)+sqrt(3.d0*hMin+2.d0*sqrt(2.d0*hMin*hMax)+sqrt(2.d0/g)*(uL-uR)*sqrt(hMin)))**2
    else
     SharpInitialGuess = sqrt(hMin*hMax)*(1.d0+(sqrt(2.d0)*(uL-uR))/(sqrt(g*hMin)+sqrt(g*hMax)))
    END if
END function SharpInitialGuess

double precision function phi(g,h,hL,hR,uL,uR)  
    !Depth function phi(h) = f(h,hL,uL) + f(h,hR,uR) + uR - uL
    !See Toro(2001) Sec 5.3
    implicit none
    double precision :: g,h,hL,hR,uL,uR
    double precision :: fL, fR
    ! Compute f(h,hL,uL)
    if (h<=hL) then 
        fL = 2.d0*(sqrt(g*h)-sqrt(g*hL)) 
    else 
        fL = (h-hL)*sqrt(0.5d0*g*(h+hL)/h/hL)
    END if
    ! Compute f(h,hR,uR) 
    if (h<=hR) then 
        fR = 2.d0*(sqrt(g*h)-sqrt(g*hR)) 
    else 
        fR = (h-hR)*sqrt(0.5d0*g*(h+hR)/h/hR)
    END if
    phi = fL + fR + uR - uL
END function phi

double precision function phip(g,h,hL,hR)
    !Derivative of depth function phi'(h) = f'(h,hL,uL) + f'(h,hR,uR)
    implicit none
    double precision :: g,h,hL,hR
    double precision :: fpL, fpR
    ! Compute f'(h,hL,uL)
    if (h<=hL) then 
        fpL = sqrt(g/h) 
    else 
        fpL = g*(2.d0*h**2+h*hL+hL**2)/(2.d0*sqrt(2.d0*g)*h**2*hL*sqrt(1/h+1/hL))
    END if
    ! Compute f'(h,hR,uR) 
    if (h<=hR) then 
        fpR = sqrt(g/h) 
    else 
        fpR = g*(2.d0*h**2+h*hR+hR**2)/(2.d0*sqrt(2.d0*g)*h**2*hR*sqrt(1/h+1/hR))
    END if
    phip = fpL + fpR  
END function phip



