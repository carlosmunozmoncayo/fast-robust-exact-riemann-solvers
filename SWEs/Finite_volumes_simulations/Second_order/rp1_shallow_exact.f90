
! =========================================================
subroutine rp1(maxmx,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =========================================================

! Solve Riemann problem for the 1D shallow water equations using 
! Ostrowski-Newton's method with two-shocks initial guess

! waves: 2
! equations: 2

! Conserved quantities:
!       1 depth (h)
!       2 momentum (hu)

    implicit real (kind=8) (a-h,o-z)

    integer, intent(in) :: maxmx, meqn, mwaves, mbc, mx, maux
    real (kind=8), dimension(meqn,1-mbc:maxmx+mbc), intent(in) :: ql, qr
    real (kind=8), dimension(maux,1-mbc:maxmx+mbc), intent(in) :: auxl, auxr
    real (kind=8), dimension(meqn, mwaves, 1-mbc:maxmx+mbc), intent(out) :: wave
    real (kind=8), dimension(meqn, 1-mbc:maxmx+mbc), intent(out) :: amdq, apdq
    real (kind=8), dimension(mwaves, 1-mbc:maxmx+mbc), intent(out) :: s

    real (kind=8) :: u_l, u_r, h_l, h_r, c_l, c_r
    real (kind=8) :: hsqrt_l, hsqrt_r, u_hat, h_hat, c_hat, grav
    real (kind=8) :: h_m, hu_m
    real (kind=8) :: h_bar, u_bar
    real (kind=8) :: cs11,cs12,cs21,cs22,ss1,ss2 !Characteristic speeds if rarefactions, shock speeds if shocks
    real (kind=8) :: rho_m, rhou_m, E_m, s1, s2
    integer :: m, i, mw
    
    integer :: switch
    !real (kind=8) :: max_wave_speed
!If switch = 0, q_bar has been found 
!	     1, both waves are rarefactions
!	     2, both waves are shocks
!	     3, 1-wave is rarefaction and 2-wave is shock
!	     4, 1-wave is shock and 2-wave is rarefaction
   
    common /cparam/  grav   
    

    tol=1.0E-12
    do i=2-mbc,mx+mbc
    switch=1 !Set to 0 if we find q_bar
        h_l = qr(1,i-1)
        h_r = ql(1,i)
        u_l = qr(2,i-1) / qr(1,i-1)
        u_r = ql(2,i  ) / ql(1,i  )
        !Left and right states sound speeds
        c_l = dsqrt(grav*h_l)
        c_r = dsqrt(grav*h_r)
        !Dry velocities
        d_l=u_l+2*c_l
        d_r=u_r-2*c_r

         !First we handle dry states
        if (h_l<tol .AND. h_r<tol) then
            wave(1,1,i) = 0
            wave(2,1,i) = 0
            wave(1,2,i) = 0
            wave(2,2,i) = 0
            s(1,i)=0
            s(2,i)=0
            amdq(1,i)=0
            amdq(2,i)=0
            apdq(1,i)=0
            apdq(2,i)=0
            cycle
        END if
    
        !We check if there is a dry state

        if (d_l<d_r) then !If a dry state is created in the middle
        !if (1.EQ.0) then
            h_m=0.
            u_m=0.
            cs11=u_l-c_l
            cs12=d_l !There is no need to make this assignation, it can be directly computed in the if
            cs21=d_r !This is also just used once in the next if 
            cs22=u_r+c_r
            s(1,i) = cs11
            s(2,i) = cs22
            if ((cs11<0) .AND. (cs12>0))   then   !If 1-wave is transonic
                u_bar=(1.0/3.0)*(u_l+2.0*c_l)
                h_bar=(u_bar*u_bar)/grav    !q_bar has been found
                switch=0
            else if ((cs21<0) .AND. (cs22>0)) then  !If 2-wave is transonic
                u_bar=(1.0/3.0)*(u_r-2.0*c_r)
                h_bar=(u_bar*u_bar)/grav   !q_bar has been found
                switch=0
            end if 
        else
            h_m=hStarWetStates(grav,h_l,h_r,u_l,u_r,tol)
            u_m=uStar(grav,h_m,h_l,h_r,u_l,u_r)
            c_m=dsqrt(grav*h_m)  
            !Determining nature of 1 wave and computing q_bar for a 1-transonic rarefaction
            !Also computing 1-speed

            if (h_m .LE. h_l) then  !If 1-wave is rarefaction
                 !We compute characteristic speeds for 1-wave
                cs11=u_l-c_l
                cs12=u_m-c_m
                s(1,i)=cs11
                if ((cs11<0) .AND. (cs12>0))   then   !If 1-wave is transonic
                    u_bar=(1.0/3.0)*(u_l+2.0*c_l)
                    h_bar=(u_bar*u_bar)/grav    !q_bar has been found
                    switch=0 
                END if
            else
                ss1= u_l-c_l*dsqrt(((h_m+h_l)*h_m)/(2*h_l*h_l))  !Shock speed for 1-wave
                s(1,i)=ss1
            END if

            !Determining nature of 2 wave and computing q_bar for 2-transonic rarefaction
            !Also computing 2-speed

            if (h_m .LE. h_r) then  !If 2-wave is rarefaction
                 !We compute characteristic speeds for 2-wave
                cs21=u_m+c_m
                cs22=u_r+c_r
                s(2,i)=cs22
                if ((cs21<0) .AND. (cs22>0))   then   !If 2-wave is transonic
                    u_bar=(1.0/3.0)*(u_r-2.0*c_r)
                    h_bar=(u_bar*u_bar)/grav    !q_bar has been found
                    switch=0
                END if
            else
                ss2= u_r+c_r*dsqrt(((h_m+h_r)*h_m)/(2*h_r*h_r))  !Shock speed for 1-wave
                s(2,i)=ss2
            END if

            !Computing q_bar if there are no transonic rarefactions
            if (switch/=0) then
                if (s(1,i)>0) then
                    h_bar=h_l
                    u_bar=u_l
                else if (s(2,i)<0) then
                    h_bar=h_r
                    u_bar=u_r
                else 
                    h_bar=h_m
                    u_bar=u_m
                END if
            END if
        END if
        
        !Computing waves
        hu_m=h_m * u_m
        wave(1,1,i) = h_m - h_l
        wave(2,1,i) = hu_m - qr(2,i-1)
        wave(1,2,i) = h_r - h_m
        wave(2,2,i) = ql(2,i) - hu_m
        
        !Computing fluctuations using flux differencing
        !----------
        amdq(1,i)=h_bar*u_bar-h_l*u_l
        amdq(2,i)=h_bar*u_bar*u_bar+0.5d0*(grav*h_bar*h_bar)-h_l*u_l*u_l-0.5d0*(grav*h_l*h_l)
        apdq(1,i)=h_r*u_r-h_bar*u_bar
        apdq(2,i)=h_r*u_r*u_r+0.5d0*(grav*h_r*h_r)-h_bar*u_bar*u_bar-0.5d0*(grav*h_bar*h_bar)
        !--------
        
    end do 

end subroutine rp1

real (kind=8) function uStar(g,hStar,hL,hR,uL,uR)
  !This is based on TORO(2001) Sec 5.3
  !Here we assume that hStar has been found 
  implicit real (kind=8) (a-h,o-z)

  ! Compute fL = f(hStar,hL,uL)
  if (hStar<=hL) then 
     fL = 2.d0*(sqrt(g*hStar)-sqrt(g*hL)) 
  else 
     fL = (hStar-hL)*sqrt(0.5d0*g*(hStar+hL)/hStar/hL)
  END if
  ! Compute fR = f(hStar,hR,uR) 
  if (hStar<=hR) then 
     fR = 2.d0*(sqrt(g*hStar)-sqrt(g*hR)) 
  else 
     fR = (hStar-hR)*sqrt(0.5d0*g*(hStar+hR)/hStar/hR)
  END if
  uStar = 0.5d0*(uL+uR)+0.5d0*(fR-fL)
END function uStar

real (kind=8) function hStarWetStates(g,hL,hR,uL,uR,tol)
    !This function consider the case where both states are wet and gives an estimate on hStar
    implicit none
    integer :: iterate, iterations_single_RP
    real (kind=8):: g,hL,hR,uL,uR,tol
    real (kind=8):: hMin, hMax, fMin, fMax
    real (kind=8):: h,phi, phip, phiR
    real (kind=8):: TwoRarefactionInitialGuess
    real (kind=8) :: yz
    real (kind=8) :: phi1, phi2, phip1

    hMin = min(hL,hR)
    fMin = phi(g,hMin,hL,hR,uL,uR)

    iterate = 1
    iterations_single_RP=0
    if (0<=fMin) then
        ! In this case both waves are rarefactions. We know the solution in this case.
        hStarWetStates = TwoRarefactionInitialGuess(g,hL,hR,uL,uR)
        iterate = 0
    else
        h=0.5d0*(hL+hR)
        h=max(hMin,(hL*yz(g,h,hL)+hR*yz(g,h,hR)-uR+uL)/(yz(g,h,hL)+yz(g,h,hR)))

        phi1=phi(g,h,hL,hR,uL,uR)
        phip1=phip(g,h,hL,hR)
        h = h-phi1/phip1
        phi2=phi(g,h,hL,hR,uL,uR)

        h=h-(phi2*phi1)/(phip1*(phi1-2*phi2))
        phiR=phi(g,h,hL,hR,uL,uR)

        iterations_single_RP=iterations_single_RP+1

        if (abs(phiR)<tol) then
            iterate=0
            hStarWetStates = h
        else
            !Extra Newton step to ensure monotone convergence
            h=max(hMin, h-phiR/phip(g,h,hL,hR))
            iterations_single_RP=iterations_single_RP+1
        end if

    END if


    if (iterate==1) then
        !We approach the solution through Newton's method
        !h_{i+1}=h_i-phi(h_i)/phi'(h_i)
        !Start iterative process
        do while(.true.)

            phiR=phi(g,h,hL,hR,uL,uR)
            if (abs(phiR) <tol) then
                exit
            end if
            h = h-phiR/phip(g,h,hL,hR)   !Updating hstar
            iterations_single_RP=iterations_single_RP+1

            if (isnan(h).or. isnan(phiR)) then
                print *, "Aborted because of Nan, hstar reached", h
                print *, "h:", h
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

real (kind=8) function yz(g, h, hz)
  implicit none
  real (kind=8) :: g, h, hz
  yz=dsqrt(0.5d0*g*(h+hz)/(h*hz))
end function yz


real (kind=8) function TwoRarefactionInitialGuess(g,hL,hR,uL,uR)
  !This the solution based on having two rarefactions
  ! Note: this can be proven to be an upper bound for hStar
  implicit real (kind=8) (a-h,o-z)
  TwoRarefactionInitialGuess = (uL-uR+2.d0*sqrt(g)*(sqrt(hL)+sqrt(hR)))**2.d0/16.d0/g

END function TwoRarefactionInitialGuess

real (kind=8) function SharpInitialGuess(g,hL,hR,uL,uR)
  !This is a sharper initial guess for hStar
  !It is designed to work for when one state is near to be dry and the other is not
  implicit real (kind=8) (a-h,o-z)

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

real (kind=8) function wMax(g,hStar,hL,hR,uL,uR)
  implicit real (kind=8) (a-h,o-z)
  w1 = uL - sqrt(g*hL)*sqrt((1.0d0 + max((hStar-hL)/2.0d0/hL,0.0d0)) * (1.0d0+max((hStar-hL)/hL,0.0d0)))
  w2 = uR + sqrt(g*hR)*sqrt((1.0d0 + max((hStar-hR)/2.0d0/hR,0.0d0)) * (1.0d0+max((hStar-hR)/hR,0.0d0)))
  wMax = max(abs(w1),abs(w2))
END function wMax

real (kind=8) function phi(g,h,hL,hR,uL,uR)  
  !Depth function phi(h) = f(h,hL,uL) + f(h,hR,uR) + uR - uL
  !See Toro(2001) Sec 5.3
  implicit real (kind=8) (a-h,o-z)
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

real (kind=8) function phip(g,h,hL,hR)
  !Derivative of depth function phi'(h) = f'(h,hL,uL) + f'(h,hR,uR)
  implicit real (kind=8) (a-h,o-z)
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




