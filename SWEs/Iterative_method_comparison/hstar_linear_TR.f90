! RS with Linear polynomials, 1 Newton step, and TR initial guess
!The linear polynomial approximates phi from below bewteen two bounds of hstar
!Its root provides an estimate of hstar from above
! =========================================================
subroutine hstar(n_data, rp_data,conv_criteria, tol, grav, exec_time)
! =========================================================

    implicit none

    integer, intent(in) :: n_data, conv_criteria
    double precision, dimension(n_data,4), intent(in) :: rp_data
    double precision, intent(in) :: tol, grav
    double precision, intent(out)  :: exec_time

    double precision :: hStarWetStates, TwoRarefactionInitialGuess
    double precision :: hStarWetStates_stagnation
    integer :: iterations,i

    double precision :: start_time, finish_time, dummy_h
 
    iterations=0

    !We use a convergence criteria based on stagnation
    call cpu_time(start_time)
    do i=1,n_data
        dummy_h=hStarWetStates(grav, rp_data(i,1), rp_data(i,2), &
                rp_data(i,3), rp_data(i,4), tol,iterations)
    end do 
    call cpu_time(finish_time)
    exec_time=finish_time-start_time
    print *,"Iterations Linear Polynomails TR: ", iterations/(n_data +0.d0)

end subroutine hstar


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

double precision function hStarWetStates(g,hL,hR,uL,uR,tol,iterations)
    !This function consider the case where both states are wet and gives an estimate on hStar 
    implicit none
    integer :: iterate, iterations, iterations_single_RP
    double precision:: g,hL,hR,uL,uR,tol
    double precision:: hMin, hMax, fMin, fMax
    double precision:: h,phi, phip, phiR
    double precision:: TwoRarefactionInitialGuess
    double precision:: hStarR, hStarL,hStarR_old,hStarL_old
    double precision:: hStarLFromQuadPhiFromAbove, hStarRFromQuadPhiFromBelow
    double precision:: hStarRFromLinearPhiFromBelow, phiL
    ! We estimate hstar from below
    hMin = min(hL,hR)
    hMax = max(hL,hR)
    fMin = phi(g,hMin,hL,hR,uL,uR)
    fMax = phi(g,hMax,hL,hR,uL,uR)


 
    ! This is meant to provide estimates of hStar from the left and right
    iterate = 1 
    

    if (0<=fMin) then 
        ! In this case both waves are rarefactions. We know the solution in this case.
        hStarWetStates = TwoRarefactionInitialGuess(g,hL,hR,uL,uR)
        iterate = 0
    else if (abs(fMax).LE. tol ) then 
        ! in case that hMax hits the solution (unlikely) 
        hStarWetStates = hMax
        iterate = 0
    else if (fMax < 0) then 
        ! In this case hMax < hStar so both waves are shocks. 
        ! We need to iterate. Here two estimates from the left and right are:
        hStarL = hMax !hMin also works but is farther away from hStar
        hStarR = TwoRarefactionInitialGuess(g,hL,hR,uL,uR) ! this estimate is guaranteed to be from the right 
    else
        ! Here we have one rarefaction and one shock 
        ! We need to iterate. Here two estimates from the left and right are: 
        hStarL = hMin
        hStarR = min(hMax,TwoRarefactionInitialGuess(g,hL,hR,uL,uR)) 
    END if

    if (iterate==1) then 
        iterations_single_RP=0
        phiR = phi(g,hStarR,hL,hR,uL,uR)
        hStarL = max(hStarL,hStarR-phiR/phip(g,hStarR,hL,hR))    
        phiL = phi(g,hStarL,hL,hR,uL,uR) 
        iterations=iterations+1
        !Start iterative process 
        do while(.true.)
        !Check if hStarL and hStarR are close enough
            if (abs(hstarL-hStarR) <tol) then
                exit
            end if
            !Save old estimates of hStar 
            !hStarL_old = hStarL
            hStarR_old = hStarR
            !Compute new estimates of hStar
            hStarR = hStarRFromLinearPhiFromBelow(g,hStarL,hStarR,hL,hR,uL,uR,phiR,phiL)
            iterations=iterations+1
            iterations_single_RP=iterations_single_RP+1
            !Evaluate depth function from the right 
            phiR = phi(g,hStarR,hL,hR,uL,uR)

            !Check if hStar or phii are NaN. 
            ! This is due to estimates from the left and from the right being the same (hStarL=hStarR)
            ! If true we are done

            if (isnan(hStarR).or.isnan(phiR)) then 
               hStarR=hStarR_old
                exit
            END if
            
            if (iterations_single_RP>40) then
                print *, "Iterations exceeded", hStarL, hStarR, hL, hR, uL, uR
                call abort
            end if
            !Check if I have reached the tolerance 
            if (abs(phiR)<tol) then 
                exit
            END if
        END do
        ! return estimation of hStar
        hStarWetStates = hStarR
    END if
END function hStarWetStates

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


double precision function hStarRFromLinearPhiFromBelow(g,hStarL,hStarR,hL,hR,uL,uR,phiR,phiL)
  ! We start considering two estimates of hStar. One from the left and one from the right. 
  ! We use these estimates to construct (a priori) a linear approximation of phi from below. 
  ! This linear approximation is monotonically increasing. 
  ! We find the root (a priori) of this linear approximation to obtain a new estimate of hStar from the right
  implicit none
  double precision :: g,hStarL,hStarR,hL,hR,uL,uR
  double precision :: Delta, phip, phi, phiDDiff2
  double precision :: phiR
  double precision :: phiL, slope

  slope = (phiR - phiL)/(hStarR - hStarL)
  hStarRFromLinearPhiFromBelow = hStarR - phiR/slope
END function hStarRFromLinearPhiFromBelow

