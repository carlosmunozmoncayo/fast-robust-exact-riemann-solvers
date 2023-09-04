! Exact RS Euler equations using quadratic interpolation with
!one Hermite polynomial and Two Rarefactions initial guess

!This module receives a matrix rp_data of Riemann problems:
!rp_data=[problem_1,problem_2,...], where
!problem_i=[pl_i,ul_i,rhol_i,pr_i,ur_i,rhor_i]
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

    !We use a convergence criteria based on the residual
    call cpu_time(start_time)
    do i=1,n_data
        a_l = dsqrt(gamma*rp_data(i,1)/rp_data(i,3))
        a_r = dsqrt(gamma*rp_data(i,4)/rp_data(i,6))
        dummy_p=pStarNonVacuum(gamma,rp_data(i,1),rp_data(i,2), &
           rp_data(i,3),rp_data(i,4),rp_data(i,5),rp_data(i,6),a_l,a_r,tol,iterations)
    END do 
    call cpu_time(finish_time)
    exec_time=finish_time-start_time
    print *,"Iterations Quadratic TR: ", iterations/(n_data +0.d0)

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
    double precision :: pMax, phiMax

    !To count iterations required for a singe Riemann Problem
    integer:: iterations_single_RP
    !To decide whether we iterate or not
    integer:: iterate


    !Our estimate for pstar
    double precision :: p, phiR, p1FromQuadPhiFromAbove, p2FromQuadPhiFromBelow
    double precision :: p1_old, p2_old, p1,p2


    !We use CL and CR instead of AL and AR to avoid confusion with sound speeds aL and aR
    CL=2/((gamma+1)*rhoL)
    CR=2/((gamma+1)*rhoR)
    BL=(gamma-1)*pL/(gamma+1)
    BR=(gamma-1)*pR/(gamma+1)


    pMin=min(pL,pR)
    pMax=max(pL,pR)
    phiMin=phi(gamma,pMin,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)
    phiMax=phi(gamma,pMax,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)

    !!!!!!!!!!!!!!!!!!!!!!!
    !!!!! ALGORITHM 1 !!!!! 
    !!!!!!!!!!!!!!!!!!!!!!!
    ! This is meant to provide estimates of pstar from the left and right
    iterate = 1 
    if (0<=phiMin) then 
        ! In this case both waves are rarefactions. We know the solution in this case.
        pStarNonVacuum = TwoRarefactionInitialGuess(gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol)
        iterate = 0
    else if (abs(phiMax).LE. tol ) then 
        ! in case that hMax hits the solution (unlikely) 
        pStarNonVacuum = pMax
        iterate = 0
    else if (phiMax < 0) then 
        ! In this case hMax < pstar so both waves are shocks. 
        ! We need to iterate. Here two estimates from the left and right are:
        p1 = pMax !hMin also works but is farther away from pstar
        p2 = TwoRarefactionInitialGuess(gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol) ! this estimate is guaranteed to be from the right 
    else
        ! Here we have one rarefaction and one shock 
        ! We need to iterate. Here two estimates from the left and right are: 
        p1 = pMin
        p2 = min(pMax,TwoRarefactionInitialGuess(gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol)) !
    END if

    if (iterate==1) then 
        ! The idea is to construct (a priori) a quadratic approx of phi from below 
        ! and find the root (a priori) of that approximation. 
        ! Each root will provide the new estimate from above of pstar, p2.
        ! We iterate on this process until some tolerance is achieved. 

        ! Before starting we improve the estimate from below via one 'classic' Newton iteration
        ! For this iteration we start with p2 which is the best estimate we have so far
        ! NOTE: due to the concavity of phi, the classic Newton iteration always estimates from below
        phiR = phi(gamma,p2,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)
        p1 = max(p1,p2-phiR/phip(gamma,p2,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR))     
        !Start iterative process 
        iterations_single_RP=0 !to know where to stop if diverges
        do while(.true.)
            !Check if p1 and p2 are close enough
            if (abs(p1-p2) <tol) then
                exit
            end if

            !Save old estimates of pstar 
            !p1_old = p1
            p2_old = p2
            !Compute new estimates of pstar
            p2 = p2FromQuadPhiFromBelow(gamma,p1,p2,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR,phiR)
            iterations=iterations+1
            iterations_single_RP=iterations_single_RP+1
            !Evaluate depth function from the right 
            phiR = phi(gamma,p2,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)
            !phiL = phi(gamma,p1,pL,uL,rhoL,pR,uR,rhoR,aL,aR)

            !Check if pstar or phii are NaN. 
            ! This is due to estimates from the left and from the right being the same (p1=p2)
            ! If true we are done
            if (isnan(p2).or.isnan(phiR)) then 
                p2=p2_old
                exit
            END if
            if (iterations_single_RP>40) then
                print *, "Iterations exceeded, pL,pR,uL,uR,rhoR,rhoL: ", pL,pR,uL,uR,rhoL,rhoR
                print *, "pstar reached: ",p2
                call abort
            END if
            !Check if I have reached the tolerance 
            if (abs(phiR)<tol) then 
                exit
            END if
        END do
        ! return estimation of pstar
        pStarNonVacuum = p2
    END if
    extra=uStar(gamma,p,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)
END function pStarNonVacuum

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

double precision function phidd112(gamma,p1,p2,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)
  ! Second divided difference phi[p1,p1,p2]
  implicit none
  double precision :: gamma,p1,p2,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR
  double precision :: phidd12, phip
  phidd112 = (phidd12(gamma,p1,p2,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR) - &
            phip(gamma,p1,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR))/(p2-p1)
END function phidd112

double precision function phidd122(gamma,p1,p2,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR) 
  ! Second divided difference phi[p1,p2,p2]
  implicit none
  double precision :: gamma,p1,p2,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR
  double precision :: phidd12, phip
  phidd122 = (phip(gamma,p2,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR) - &
            phidd12(gamma,p1,p2,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR))/(p2-p1)
END function phidd122

double precision function phidd12(gamma,p1,p2,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)
    ! First divided difference phi[p1,p2]
    implicit none
    double precision :: gamma,p1,p2,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR
    double precision :: phip, phi
    if (p1 == p2) then
        phidd12= phip(gamma,p1,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)
    else
        phidd12 = (phi(gamma,p2,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR) - &
                phi(gamma,p1,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR))/(p2-p1)
    end if
END function phidd12


double precision function p2FromQuadPhiFromBelow(gamma,p1,p2,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR,phiR)
  ! We start considering two estimates of pStar. One from the left(p1) and one from the right(p2). 
  ! We use these estimates to construct (a priori) a quadratic approximation of phi from below. 
  ! This quad approximation is monotonically increasing and concave down (just as phi). 
  ! We find the root (a priori) of this quadratic approximation to obtain a new estimate of pStar from the right
  implicit none
  double precision :: gamma,p1,p2,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR
  double precision :: Delta, phip, phi , phidd122
  double precision :: phiR
  
  Delta = phip(gamma,p2,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)**2.d0 - &
          4.d0*phiR*phidd122(gamma,p1,p2,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)
  if (Delta<0) then 
     print *, 'Delta < 0 when computing the root of the quad approx of phi from above'
     print *, "pL: ", pL, " pR: ", pR
     print *, "uL: ", uL, " uR: ", uR
     print *, 'Delta: ', Delta
     call abort
  END if
  p2FromQuadPhiFromBelow = p2 - 2.d0*phiR/  &
                           (phip(gamma,p2,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)+dsqrt(Delta))

END function p2FromQuadPhiFromBelow
