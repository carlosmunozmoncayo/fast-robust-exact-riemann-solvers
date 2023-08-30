! RS with Quadratic polynomials, 1 Newton step, and TR initial guess
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
    print *,"Iterations Quadratic Polynomails TR: ", iterations/(n_data +0.d0)

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
        hStarR = min(hMax,TwoRarefactionInitialGuess(g,hL,hR,uL,uR)) !hMax is enough but let's try to get closer to hStar
    END if

    if (iterate==1) then 
        iterations_single_RP=0
        ! The idea is to construct (a priori) a quadratic approx of phi from above and from below of phi
        ! and find the root (a priori) of those approximations. 
        ! Each root will provide new estimates of hStarL and hStarR.
        ! We iterate on this process until some tolerance is achieved. 
        ! Finally, we use the estimate from the right to assure positivity. 

        ! Before starting we improve the estimate from below via one 'classic' Newton iteration
        ! For this iteration we start with hStarR which is the best estimate we have so far
        ! NOTE: due to the concavity of phi, the classic Newton iteration always estimates from below
        phiR = phi(g,hStarR,hL,hR,uL,uR)
        hStarL = max(hStarL,hStarR-phiR/phip(g,hStarR,hL,hR))     
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
            hStarR = hStarRFromQuadPhiFromBelow(g,hStarL,hStarR_old,hL,hR,uL,uR,phiR)
            iterations=iterations+1
            iterations_single_RP=iterations_single_RP+1
            !Evaluate depth function from the right 
            phiR = phi(g,hStarR,hL,hR,uL,uR)
            !phiL = phi(g,hStarR,hL,hR,uL,uR)

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

double precision function phiDDiff1(g,hStarL,hStarR,hL,hR,uL,uR)
  ! Second divided difference phi[hStarL,hStarL,hStarR]
  implicit none
  double precision :: g,hStarL,hStarR,hL,hR,uL,uR
  double precision :: phiDiff,phip
  phiDDiff1 = (phiDiff(g,hStarL,hStarR,hL,hR,uL,uR) - phip(g,hStarL,hL,hR))/(hStarR-hStarL)
END function phiDDiff1

double precision function phiDDiff2(g,hStarL,hStarR,hL,hR,uL,uR) 
  ! Second divided difference phi[hStarL,hStarR,hStarR]
  implicit none
  double precision :: g,hStarL,hStarR,hL,hR,uL,uR
  double precision :: phip, phiDiff
  phiDDiff2 = (phip(g,hStarR,hL,hR) - phiDiff(g,hStarL,hStarR,hL,hR,uL,uR))/(hStarR-hStarL)
END function phiDDiff2

double precision function phiDiff(g,hStarL,hStarR,hL,hR,uL,uR)
  ! First divided difference phi[hStarL,hStarR]
  implicit none
  double precision :: g,hStarL,hStarR,hL,hR,uL,uR
  double precision :: phip, phi
  if (hstarL == hstarR) then
   phiDiff= phip(g,hstarL,hL,hR)
  else
  phiDiff = (phi(g,hStarR,hL,hR,uL,uR) - phi(g,hStarL,hL,hR,uL,uR))/(hStarR-hStarL)
  end if
END function phiDiff

double precision function hStarLFromQuadPhiFromAbove(g,hStarL,hStarR,hL,hR,uL,uR)
  ! We start considering two estimates of hStar. One from the left and one from the right. 
  ! We use these estimates to construct (a priori) a quadratic approximation of phi from above. 
  ! This quad approximation is monotonically increasing and concave down (just as phi). 
  ! We find the root (a priori) of this quadratic approximation to obtain a new estimate of hStar from the left
  implicit none
  double precision :: g,hStarL,hStarR,hL,hR,uL,uR
  double precision :: Delta, phip, phi, phiDDiff1

  Delta = phip(g,hStarL,hL,hR)**2.d0 - 4.d0*phi(g,hStarL,hL,hR,uL,uR)*phiDDiff1(g,hStarL,hStarR,hL,hR,uL,uR)
  if (Delta<0) then 
     print *, "hstarL_old: ", hStarL, " hstarR_old: ", hStarR
     print *, "hL: ", hL, " hR: ", hR
     print *, "uL: ", uL, " uR: ", uR
     print *, 'Delta < 0 when computing the root of the quad approx of phi from above'
     print *, 'Delta: ', Delta
     call abort
  END if
  hStarLFromQuadPhiFromAbove = hStarL - 2.d0*phi(g,hStarL,hL,hR,uL,uR)/(phip(g,hStarL,hL,hR)+sqrt(Delta))
END function hStarLFromQuadPhiFromAbove

double precision function hStarRFromQuadPhiFromBelow(g,hStarL,hStarR,hL,hR,uL,uR,phiR)
  ! We start considering two estimates of hStar. One from the left and one from the right. 
  ! We use these estimates to construct (a priori) a quadratic approximation of phi from below. 
  ! This quad approximation is monotonically increasing and concave down (just as phi). 
  ! We find the root (a priori) of this quadratic approximation to obtain a new estimate of hStar from the right
  implicit none
  double precision :: g,hStarL,hStarR,hL,hR,uL,uR
  double precision :: Delta, phip, phi, phiDDiff2
  double precision :: phiR
  
  !Delta = phip(g,hStarR,hL,hR)**2.d0 - 4.d0*phi(g,hStarR,hL,hR,uL,uR)*phiDDiff2(g,hStarL,hStarR,hL,hR,uL,uR)
  Delta = phip(g,hStarR,hL,hR)**2.d0 - 4.d0*phiR*phiDDiff2(g,hStarL,hStarR,hL,hR,uL,uR)
  if (Delta<0) then 
     print *, "hstarL_old: ", hStarL, " hstarR_old: ", hStarR
     print *, "hL: ", hL, " hR: ", hR
     print *, "uL: ", uL, " uR: ", uR
     print *, 'Delta < 0 when computing the root of the quad approx of phi from below'
     print *, 'Delta: ', Delta
     call abort
  END if
  hStarRFromQuadPhiFromBelow = hStarR - 2.d0*phiR/(phip(g,hStarR,hL,hR)+sqrt(Delta))

END function hStarRFromQuadPhiFromBelow
