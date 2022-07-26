!RS using Newton's method with QA initial guess
! =========================================================
subroutine hstar(n_data, rp_data,conv_criteria, tol, grav, exec_time)
! =========================================================

    implicit none

    integer, intent(in) :: n_data, conv_criteria
    double precision, dimension(n_data,4), intent(in) :: rp_data
    double precision, intent(in) :: tol, grav
    double precision, intent(out)  :: exec_time

    double precision :: hStarWetStates, TwoRarefactionInitialGuess
    double precision :: SharpInitialGuess
    double precision :: hStarWetStates_stagnation
    integer :: iterations,i

    double precision :: start_time, finish_time, dummy_h
 
    iterations=0
    if (conv_criteria==1) then
        !We use a convergence criteria based on the residual
        call cpu_time(start_time)
        do i=1,n_data
            dummy_h=hStarWetStates(grav, rp_data(i,1), rp_data(i,2), &
                    rp_data(i,3), rp_data(i,4), tol,iterations)
        end do 
        call cpu_time(finish_time)
        exec_time=finish_time-start_time
        print *,"Iterations Newton QA: ", iterations/(n_data +0.d0)

    else if (conv_criteria==0 ) then
        !We use a convergence criteria based on the residual
        call cpu_time(start_time)
        do i=1,n_data
            dummy_h=hStarWetStates_stagnation(grav, rp_data(i,1), rp_data(i,2), &
                    rp_data(i,3), rp_data(i,4), tol,iterations)
        end do 
        call cpu_time(finish_time)
        exec_time=finish_time-start_time
        print *,"Iterations Newton QA: ", iterations/(n_data +0.d0)
    else
        print *, "Convergence criteria not properly specified"
        call abort

    END if
end subroutine hstar

double precision function hStarWetStates(g,hL,hR,uL,uR,tol,iterations)
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

    iterate = 1 
    iterations_single_RP=0

    if (0<=fMin) then 
        ! In this case both waves are rarefactions. We know the solution in this case.
        hStarWetStates = TwoRarefactionInitialGuess(g,hL,hR,uL,uR)
        iterate = 0
    else
        h=SharpInitialGuess(g,hL,hR,uL,uR)
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
                print *, "Aborted because of Nan, hstar reached", h
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

double precision function hStarWetStates_stagnation(g,hL,hR,uL,uR,tol,iterations)
    !This function consider the case where both states are wet and gives an estimate on hStar 
    implicit none
    integer :: iterate, iterations, iterations_single_RP
    double precision:: g,hL,hR,uL,uR,tol
    double precision:: hMin, hMax, fMin, fMax
    double precision:: h, hold,phi, phip, phiR
    double precision:: TwoRarefactionInitialGuess, SharpInitialGuess

    ! We estimate hstar from below
    hMin = min(hL,hR)
    ! hMax = max(hL,hR)
    fMin = phi(g,hMin,hL,hR,uL,uR)
    ! fMax = phi(g,hMax,hL,hR,uL,uR)

    iterate = 1 
    iterations_single_RP=0

    
    if (0<=fMin) then 
        ! In this case both waves are rarefactions. We know the solution in this case.
        hStarWetStates_stagnation = TwoRarefactionInitialGuess(g,hL,hR,uL,uR)
        iterate = 0
    else
        h=SharpInitialGuess(g,hL,hR,uL,uR)
        hold=h
    END if

    if (iterate==1) then 
        !We approach the solution through Newton's method
        !h_{i+1}=h_i-phi(h_i)/phi'(h_i)
        h=max(hMin, h-phi(g,h,hL,hR,uL,uR)/phip(g,h,hL,hR))
        iterations_single_RP=iterations_single_RP+1
        iterations=iterations+1
        !Start iterative process 
        do while(.true.)
            !Trying to fix 
            if (abs(2.d0*(hold-h)/(hold+h)) <tol) then
                exit
            end if
            hold=h
            h = h-phi(g,h,hL,hR,uL,uR)/phip(g,h,hL,hR)   !Updating hstar
            iterations=iterations+1
            iterations_single_RP=iterations_single_RP+1
            if (isnan(h).or. isnan(phiR)) then
                print *, "Aborted because of Nan, hstar reached", h
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
        hStarWetStates_stagnation = h
    END if  
END function hStarWetStates_stagnation

double precision function uStar(g,h,hL,hR,uL,uR)
    !This is based on TORO(2001) Sec 5.3
    !Here we assume that hStar has been found 
    implicit none
    double precision :: g,h,hL,hR,uL,uR
    double precision :: fL, fR

    ! Compute fL = f(hStar,hL,uL)
    if (h<=hL) then 
         fL = 2.d0*(sqrt(g*h)-sqrt(g*hL)) 
    else 
         fL = (h-hL)*sqrt(0.5d0*g*(h+hL)/h/hL)
    END if
    ! Compute fR = f(hStar,hR,uR) 
    if (h<=hR) then 
         fR = 2.d0*(sqrt(g*h)-sqrt(g*hR)) 
    else 
         fR = (h-hR)*sqrt(0.5d0*g*(h+hR)/h/hR)
    END if
    uStar = 0.5d0*(uL+uR)+0.5d0*(fR-fL)
END function uStar

double precision function TwoRarefactionInitialGuess(g,hL,hR,uL,uR)
    !This the solution based on having two rarefactions
    ! Note: this can be proven to be an upper bound for hStar
    implicit none
    double precision :: g,hL,hR,uL,uR
    TwoRarefactionInitialGuess = (uL-uR+2.d0*sqrt(g)*(sqrt(hL)+sqrt(hR)))**2.d0/16.d0/g
END function TwoRarefactionInitialGuess

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
