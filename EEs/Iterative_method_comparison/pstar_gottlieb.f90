! Exact RS for Euler equations using Gottlieb & Groth's method, iterating with ustar
!The stopping criteria is stagnation

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

    double precision :: pStarNonVacuum
    double precision :: a_l, a_r
    integer :: iterations, i

    double precision :: start_time, finish_time, dummy_p

    iterations=0
    
    call cpu_time(start_time)
    !This loop solves each Riemann problem in rp_data, 
    !The solutions are not stored to measure execution time
    do i=1,n_data
            a_l = dsqrt(gamma*rp_data(i,1)/rp_data(i,3))
            a_r = dsqrt(gamma*rp_data(i,4)/rp_data(i,6))
            dummy_p=pStarNonVacuum(gamma,rp_data(i,1),rp_data(i,2), &
                        rp_data(i,3),rp_data(i,4),rp_data(i,5),rp_data(i,6),a_l,a_r,tol,iterations)
            !print *, "Solution:", dummy_p

    end do 
    call cpu_time(finish_time)
    exec_time=finish_time-start_time
    
    print *,"Iterations Gottlieb:", (iterations/(0.d0+n_data))

end subroutine pstar

double precision function TwoRarefactionInitialGuess(gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol) 
  !This the exact solution for u_star based on having two rarefactions
  implicit none
  double precision :: gamma, pL, uL, rhoL, pR, uR, rhoR, aL, aR, tol
  double precision :: z, uL_t, uR_t 
  z=(aR/aL)*(pL/pR)**(0.5d0*(gamma-1)/gamma)
  uL_t=uL+2.d0*aL/(gamma-1)
  uR_t=uR-2.d0*aR/(gamma-1)
  TwoRarefactionInitialGuess = (uL_t*z+uR_t)/(1+z)

END function TwoRarefactionInitialGuess

double precision function TwoRarefactionInitialGuess_pressure(gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol) 
    !This the exact solution for p_star based on having two rarefactions
    implicit none
    double precision, intent (in) :: gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol
    double precision :: gamma1
    gamma1=gamma-1.d0
    TwoRarefactionInitialGuess_pressure = ((aL+aR-0.5d0*gamma1*(uR-uL))/(aL/pL**(gamma1/(2*gamma)) &
        +aR/pR**(gamma1/(2*gamma))))**(2*gamma/gamma1)

END function TwoRarefactionInitialGuess_pressure

double precision function pStarNonVacuum(gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol,iterations) 
    !This function estimates pstar exactly if there are two rarefactions and 
    !approximately using Gottlieb's method otherwise
    
    implicit none
    double precision :: gamma, pL, uL, rhoL, pR, uR, rhoR, aL, aR, tol
    integer:: iterations, iterations_single_RP , iterate

    double precision :: VL, VR, BL, BR

    double precision :: phi, TwoRarefactionInitialGuess

    double precision :: TwoRarefactionInitialGuess_pressure, pMin, phiMin

    double precision :: pstarL, pstarR, uStar, pstarR_p, pstarL_p, astarL, astarR
    double precision :: WL, WR, CL, CR

    iterate = 1 
    iterations_single_RP= 0

    CL=gamma*pL/aL
    CR=gamma*pR/aR

    !We use VL and VR instead of AL and AR to avoid confusion with sound speeds aL and aR
    VL=2/((gamma+1)*rhoL)
    VR=2/((gamma+1)*rhoR)

    BL=(gamma-1)*pL/(gamma+1)
    BR=(gamma-1)*pR/(gamma+1)

    !We solve the problem exactly or use the Two Rarefaction initial guess propossed by Gottlieb 

    ustar=TwoRarefactionInitialGuess(gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol)

    !Since we are working with gamma_l=gamma_r, this is the exact solution if we have two rarefactions
    !To have a fair comparison in performance with other solvers, we rule out this scenario first. If 
    !there are two rarefactions, we just compute pstar exactly

    pMin=min(pL,pR)
    phiMin=phi(gamma,pMin,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)

    if (phiMin>0) then
        pStarNonVacuum=TwoRarefactionInitialGuess_pressure(gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol)
        iterate=0
    END if   

    if (iterate==1) then 

        !This solver is crashing with this    9.4103942987765286E-003  -62.246899203685381       0.68541728873344243        
        !93.182195202653389        3.7580655775098395E-003  0.85831123762303996     

        if(ustar .LE. uL) then
            WL=0.25d0*(gamma+1)*(ustar-uL)/aL-(1.d0+(0.25d0*(gamma+1)*(ustar-uL)/aL)**2.d0)**0.5d0
            pstarL=pL+CL*(ustar-uL)*WL
            pstarL_p=(2.d0*CL*WL**3.d0)/(1.d0+WL**2.d0)
           !astarL=aL*(((gamma+1.d0)+(gamma-1.d0)*pstarL/pL)/((gamma+1)+(gamma-1)*pL/pstarL))**0.5d0
        else
            astarL=aL-0.5d0*(gamma-1.d0)*(ustar-uL)
            pstarL=pL*(astarL/aL)**(2.d0*gamma/(gamma-1.d0))
            pstarL_p=-gamma*pstarL/astarL
            
        end if

        if (ustar .GE. uR) then
            WR=0.25d0*(gamma+1)*(ustar-uR)/aR+(1.d0+(0.25d0*(gamma+1.d0)*(ustar-uR)/aR)**2.d0)**0.5d0
            pstarR=pR+CR*(ustar-uR)*WR
            pstarR_p=(2.d0*CR*WR**3.d0)/(1.d0+WR**2.d0)
            !astarR=aR*(((gamma+1.d0)+(gamma-1.d0)*pstarR/pR)/((gamma+1.d0)+(gamma-1.d0)*pR/pstarR))**0.5d0
        else
            astarR=aR+0.5d0*(gamma-1.d0)*(ustar-uR)
            pstarR=pR*(astarR/aR)**(2.d0*gamma/(gamma-1.d0))
            pstarR_p=gamma*pstarR/astarR
        end if


        ustar=ustar-(pstarL-pstarR)/(pstarL_p-pstarR_p)
        iterations=iterations+1
         
        do while(.true.)

            !Stagnation criteria proposed by Gottlieb |pstarR-pstarL|/pstarR
            if (abs(1-pstarL/pstarR)<tol) then
                exit
            END if

            !Residual based convergence criteria|phi(pstar)|< tol 
            !Would be too inefficient since requires an additional evaluation of phi

            !phiR = phi(gamma,p_star,pL,uL,rhoL,pR,uR,rhoR,aL,aR,VL,VR,BL,BR)
            !if (abs(phiR)<tol) then
            !   exit
            !END if

            if(ustar .LE. uL) then
                WL=0.25d0*(gamma+1)*(ustar-uL)/aL-(1.d0+(0.25d0*(gamma+1)*(ustar-uL)/aL)**2.d0)**0.5d0
                pstarL=pL+CL*(ustar-uL)*WL
                pstarL_p=(2.d0*CL*WL**3.d0)/(1.d0+WL**2.d0)
                !astarL=aL*(((gamma+1.d0)+(gamma-1.d0)*pstarL/pL)/((gamma+1)+(gamma-1)*pL/pstarL))**0.5d0
            else
                astarL=aL-0.5d0*(gamma-1.d0)*(ustar-uL)
                pstarL=pL*(astarL/aL)**(2.d0*gamma/(gamma-1.d0))
                pstarL_p=-gamma*pstarL/astarL
            
            end if

            if (ustar .GE. uR) then
                WR=0.25d0*(gamma+1)*(ustar-uR)/aR+(1.d0+(0.25d0*(gamma+1.d0)*(ustar-uR)/aR)**2.d0)**0.5d0
                pstarR=pR+CR*(ustar-uR)*WR
                pstarR_p=(2.d0*CR*WR**3.d0)/(1.d0+WR**2.d0)
                !astarR=aR*(((gamma+1.d0)+(gamma-1.d0)*pstarR/pR)/((gamma+1.d0)+(gamma-1.d0)*pR/pstarR))**0.5d0
            else
                astarR=aR+0.5d0*(gamma-1.d0)*(ustar-uR)
                pstarR=pR*(astarR/aR)**(2.d0*gamma/(gamma-1.d0))
                pstarR_p=gamma*pstarR/astarR
            end if

            ustar=ustar-(pstarL-pstarR)/(pstarL_p-pstarR_p)
            iterations=iterations+1
            iterations_single_RP=iterations_single_RP+1

            if (isnan(pstarL).or.isnan(ustar)) then 
                print *, "aborted at", iterations_single_RP, "iterations, nan p_star or phiR"
                print *, pL,uL,rhoL,pR,uR,rhoR, ustar, pstarL
                call abort
            END if        
        
            if (iterations_single_RP>40) then
                print *, "aborted at", iterations_single_RP, "iterations"
                print *, pL,uL,rhoL,pR,uR,rhoR
                call abort
            
            end if
        END do
     ! return estimation of pStar
     pStarNonVacuum = pstarL
    END if
END function pStarNonVacuum

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
