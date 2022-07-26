! Exact RS Euler Equations using Newton's method with several initial guess
! =========================================================

module Initial_Guesses
   implicit none
   contains

   double precision function TwoRarefactionInitialGuess&
   (gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol,CL,CR,BL,BR, pMin)
     !This the exact solution based on having two rarefactions
     double precision, intent(in) :: gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol,CL,CR,BL,BR,pMin
     double precision :: gamma1
     gamma1=gamma-1.d0
     TwoRarefactionInitialGuess = ((aL+aR-0.5d0*gamma1*(uR-uL))/(aL/pL**(gamma1/(2*gamma))  &
                                  +aR/pR**(gamma1/(2*gamma))))**(2*gamma/gamma1)
     return

   END function TwoRarefactionInitialGuess

   double precision function gs(gamma, p, CS, BS)
      double precision :: gamma, p, CS, BS
      gs=dsqrt(CS/(p+BS))
   end function gs

   double precision function TwoShockInitialGuess&
   (gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol,CL,CR,BL,BR, pMin) 
      double precision, intent(in) :: gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol,CL,CR,BL,BR,pMin
      double precision :: gs, ppv

      ppv=max(pMin,0.5d0*(pL+pR)-0.125d0*(uR-uL)*(rhoL+rhoR)*(aL+aR))
      TwoShockInitialGuess = (gs(gamma, ppv, CL, BL)*pL+gs(gamma, ppv, CR, BR)*pR-(uR-uL))/&
             (gs(gamma,ppv,CL,BL)+gs(gamma,ppv,CR,BR))
      return

   END function TwoShockInitialGuess

   double precision function AverageInitialGuess&
   (gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol,CL,CR,BL,BR, pMin)
      double precision, intent(in) :: gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol,CL,CR,BL,BR,pMin
      AverageInitialGuess=0.5d0*(pL+pR)

      return
   END function

   double precision function LinearizedSolver&
   (gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol,CL,CR,BL,BR, pMin)
      double precision, intent(in) :: gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol,CL,CR,BL,BR, pMin
      LinearizedSolver=max(pMin,0.5d0*(pL+pR)-0.125d0*(uR-uL)*(rhoL+rhoR)*(aL+aR))

      return
   END function



end module Initial_Guesses


subroutine pstar(n_data, rp_data, conv_criteria, initial_guess, tol ,exec_time)
! =========================================================
   use Initial_Guesses
   
   implicit none

   integer, intent(in) :: n_data
   double precision, dimension(n_data,6), intent(in) :: rp_data
   integer, intent(in) :: conv_criteria !If 1 we determine converge with residual, if 0 we do it with stagnation
   integer, intent(in) :: initial_guess !If 0 TR, if 1 TS, 2 if AV, 3 if PPV
   double precision, intent(in) :: tol
   double precision, intent(out) :: exec_time
    !double precision , intent(out):: salida
    
   double precision :: uStar, pStarNonVacuum!, TwoRarefactionInitialGuess, SharpInitialGuess
   double precision :: phi, phip, phidd112, phidd122, phidd12, wMax
   double precision :: u_l, u_r, p_l, p_r,rho_l,rho_r, E_l, E_r, a_l, a_r, d_l, d_r
   double precision :: hsqrt_l, hsqrt_r, u_hat, h_hat, c_hat, grav, gamma, gamma1
   double precision :: u_m, p_m, hu_m, c_m, p_bar
   double precision :: h_bar, u_bar
   double precision :: cs11,cs12,cs31,cs32,ss1,ss3 !Characteristic speeds if rarefactions, shock speeds if shocks
   double precision :: rho_m, rho_bar, rho_ml, rho_mr, rhou_m, E_m, s1, s2, S_l

   !double precision :: TwoRarefactionInitialGuess, TwoShockInitialGuess, LinearizedSolver, AverageInitialGuess, gs


   logical :: vacuum_middle
   integer :: m, i, mw, iterations

   integer :: switch

   double precision :: start_time, finish_time, dummy_p

   abstract interface
      function funcIG (gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol,CL,CR,BL,BR, pMin)
         double precision:: funcIG
         double precision, intent (in) :: gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol,CL,CR,BL,BR, pMin
      end function funcIG
   end interface
   

   procedure (funcIG), pointer :: Initial_Guess_ptr => null ()

   gamma=1.4
   gamma1=gamma-1.d0

   iterations=0


   
   if (initial_guess==0) then
      Initial_Guess_ptr => TwoRarefactionInitialGuess
   else if (initial_guess==1) then
      Initial_Guess_ptr => TwoShockInitialGuess
   else if(initial_guess==2) then
      Initial_Guess_ptr => AverageInitialGuess
   else
      Initial_Guess_ptr => LinearizedSolver
   END if

   call cpu_time(start_time)
   do i=1,n_data
         a_l = dsqrt(gamma*rp_data(i,1)/rp_data(i,3))
         a_r = dsqrt(gamma*rp_data(i,4)/rp_data(i,6))
         dummy_p=pStarNonVacuum(gamma,rp_data(i,1),rp_data(i,2), &
                     rp_data(i,3),rp_data(i,4),rp_data(i,5),rp_data(i,6),a_l,a_r,tol,iterations)

   end do 
   call cpu_time(finish_time)
   exec_time=finish_time-start_time

   print *,"Iterations Newton TS: ", iterations/(0.d0+n_data)

end subroutine pstar



double precision function pStarNonVacuum(gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol,iterations) 
   !This function consider the case where pressure is positive in the initial states and its positivity is preserved  in the middle and gives an estimate on hStar 
   implicit none
   double precision :: gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol,CL,CR,BL,BR,pMin
   double precision :: fMin, p2, p2old, phiR, extra, phi, phip, ustar
   integer :: iterate
   double precision :: Initial_Guess_ptr, TwoRarefactionInitialGuess
   double precision :: AverageInitialGuess, TwoShockInitialGuess, LinearizedSolver
   double precision :: gs

   integer:: iterations, iterations2

    !We use CL and CR instead of AL and AR to avoid confusion with sound speeds aL and aR
   CL=2/((gamma+1)*rhoL)
   CR=2/((gamma+1)*rhoR)
   BL=(gamma-1)*pL/(gamma+1)
   BR=(gamma-1)*pR/(gamma+1)

   ! We need to start with two states of h. One to the left and one to the right of hStar
   pMin = min(pL,pR)
   !pMax = max(pL,pR)
   fMin = phi(gamma,pMin,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)
   !fMax = phi(gamma,pMax,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)

   ! This is meant to provide estimates of hStar from the left and right
   iterate = 1 
   iterations2= 0
   if (0<=fMin) then 
     ! In this case both waves are rarefactions. We know the solution in this case.
     pStarNonVacuum = TwoRarefactionInitialGuess &
     (gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol,CL,CR,BL,BR, pMin)
     iterate = 0
   else
     p2=Initial_Guess_ptr(gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol,CL,CR,BL,BR, pMin)

   END if

   if (iterate==1) then 
     ! Before starting we improve the estimate from below via one 'classic' Newton iteration
     ! For this iteration we start with p2 which is the best estimate we have so far
     ! NOTE: due to the concavity of phi, the classic Newton iteration always estimates from below
     p2 = max(pMin,p2-phi(gamma,p2,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)/&
         phip(gamma,p2,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR))     
     iterations=iterations+1
     !Start iterative process 
     do while(.true.)

        !Check if p1 and p2 are close enough
        !if (abs(p1-p2) <tol) then
        !   exit
        !end if

        phiR = phi(gamma,p2,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)
        if (abs(phiR)<tol) then 
           exit
        END if
        p2old = p2

        !Compute new estimates of pStar

        p2 = p2-phiR/&
            phip(gamma,p2,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)
        iterations=iterations+1
        iterations2=iterations2+1

        if (isnan(p2).or.isnan(phiR)) then 
            print *, "aborted at", iterations2, "iterations by Nan"
            print *, pL,uL,rhoL,pR,uR,rhoR
            print*, "pstar reached: ", p2
            print*, "phiR, phipR:", phiR, phip(gamma,p2,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)
            call abort
        END if
        
        if (iterations2>40) then
            print *, "aborted at", iterations2, "iterations"
            print *, pL,uL,rhoL,pR,uR,rhoR
            print*, "pstar reached: ", p2
            print*, "phiR, phipR:", phiR, phip(gamma,p2,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)
            call abort
            
        end if
        
     END do
     ! return estimation of pStar
     pStarNonVacuum = p2
   END if
      extra=uStar(gamma,p2,pL,uL,rhoL,pR,uR,rhoR,aL,aR)
END function pStarNonVacuum


double precision function uStar(gamma,pStar,pL,uL,rhoL,pR,uR,rhoR,aL,aR)
  !This is based on TORO(2008) Sec 4.2, Prop 4.1
  !Here we assume that pStar has been found 
  implicit double precision (a-h,o-z)

  ! Compute fL = f(hStar,hL,uL)
  CL=2/((gamma+1)*rhoL)
  CR=2/((gamma+1)*rhoR)
  BL=(gamma-1)*pL/(gamma+1)
  BR=(gamma-1)*pR/(gamma+1)
  ! Compute fL(p,WL)
  if (pStar<=pL) then !Rarefaction 
     fL = (2*aL/(gamma-1))*((pStar/pL)**((gamma-1)/(2*gamma))-1)
  else 
     fL = (pStar-pL)*dsqrt(CL/(pStar+BL))
  END if
  ! Compute fR(p,WR) 
  if (pStar<=pR) then 
     fR = (2*aR/(gamma-1))*((pStar/pR)**((gamma-1)/(2*gamma))-1)
  else 
     fR = (pStar-pR)*dsqrt(CR/(pStar+BR))
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


! !Initial Guesses
! double precision function TwoRarefactionInitialGuess(gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol,CL,CR,BL,BR, pMin)
!   !This the exact solution based on having two rarefactions
!   double precision, intent(in) :: gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol,CL,CR,BL,BR,pMin
!   double precision :: gamma1
!   gamma1=gamma-1.d0
!   TwoRarefactionInitialGuess = ((aL+aR-0.5d0*gamma1*(uR-uL))/(aL/pL**(gamma1/(2*gamma))  &
!                                +aR/pR**(gamma1/(2*gamma))))**(2*gamma/gamma1)
!   return
! END function TwoRarefactionInitialGuess



! double precision function TwoShockInitialGuess(gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol,CL,CR,BL,BR, pMin) 
!    double precision, intent(in) :: gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol,CL,CR,BL,BR,pMin
!    double precision :: gs, ppv

!    ppv=max(pMin,0.5d0*(pL+pR)-0.125d0*(uR-uL)*(rhoL+rhoR)*(aL+aR))
!    TwoShockInitialGuess = (gs(gamma, ppv, CL, BL)*pL+gs(gamma, ppv, CR, BR)*pR-(uR-uL))/&
!           (gs(gamma,ppv,CL,BL)+gs(gamma,ppv,CR,BR))
!    return
! END function TwoShockInitialGuess


! double precision function gs(gamma, p, CS, BS)
!    double precision :: gamma, p, CS, BS
!    gs=dsqrt(CS/(p+BS))
!    return
! end function gs


! double precision function AverageInitialGuess(gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol,CL,CR,BL,BR, pMin)
!    double precision, intent(in) :: gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol,CL,CR,BL,BR,pMin
!    AverageInitialGuess=0.5d0*(pL+pR)
!    return
! END function



! double precision function LinearizedSolver(gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol,CL,CR,BL,BR, pMin)
!    double precision, intent(in) :: gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol,CL,CR,BL,BR, pMin
!    LinearizedSolver=max(pMin,0.5d0*(pL+pR)-0.125d0*(uR-uL)*(rhoL+rhoR)*(aL+aR))
!    return
! END function
