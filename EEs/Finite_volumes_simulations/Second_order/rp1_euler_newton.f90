
! =========================================================
subroutine rp1(maxmx,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =========================================================

! Solve Riemann problem for the 1D Euler equations
! using positive Newton's method with two-shock initial guess

! waves: 2
! equations:3

! Conserved quantities:
!       1 depth (h)
!       2 momentum (hu)

    implicit none

    integer, intent(in) :: maxmx, meqn, mwaves, mbc, mx, maux
    double precision, dimension(meqn,1-mbc:maxmx+mbc), intent(in) :: ql, qr
    double precision, dimension(maux,1-mbc:maxmx+mbc), intent(in) :: auxl, auxr
    double precision, dimension(meqn, mwaves, 1-mbc:maxmx+mbc), intent(out) :: wave
    double precision, dimension(meqn, 1-mbc:maxmx+mbc), intent(out) :: amdq, apdq
    double precision, dimension(mwaves, 1-mbc:maxmx+mbc), intent(out) :: s
    
    double precision :: uStar, pStarNonVacuum, TwoRarefactionInitialGuess, SharpInitialGuess
    double precision :: phi, phip, phidd112, phidd122, phidd12, wMax, p1FromQuadPhiFromAbove, p2FromQuadPhiFromBelow

    double precision :: u_l, u_r, p_l, p_r,rho_l,rho_r, E_l, E_r, a_l, a_r, d_l, d_r
    double precision :: hsqrt_l, hsqrt_r, u_hat, h_hat, c_hat, grav, gamma, gamma1
    double precision :: u_m, p_m, hu_m, c_m, p_bar
    double precision :: h_bar, u_bar
    double precision :: cs11,cs12,cs31,cs32,ss1,ss3 !Characteristic speeds if rarefactions, shock speeds if shocks
    double precision :: rho_m, rho_bar, rho_ml, rho_mr, rhou_m, E_m, s1, s2, S_l
    double precision :: tol

    logical :: vacuum_middle
    integer :: m, i, mw, iterations
    
    integer :: switch
    !double precision :: max_wave_speed
!If switch < 5 there are not transonic rarefactions, q bar has not been found
!	     1, both waves are rarefactions
!	     2, both waves are shocks
!	     3, 1-wave is rarefaction and 2-wave is shock
!	     4, 1-wave is shock and 2-wave is rarefaction
!If switch >= 5 there is a transonic rarefaction, q bar has been found
!            5, 1-transonic, 2-rarefaction
!            6, 1-rarefaction, 2-transonic
!            7, 1-transonic, 2-shock
!            8, 1-shock, 2-transonic
    		
    common /cparam/   gamma, gamma1


    tol=1.0E-12
    
    do i=2-mbc,mx+mbc
        !Density
        rho_l=qr(1,i-1)
        rho_r=ql(1,i  )
        ! Velocity
        u_l = qr(2,i-1) / qr(1,i-1)
        u_r = ql(2,i  ) / ql(1,i  )
        ! Pressure
        p_l = gamma1 * (qr(3,i-1) - rho_l*u_l**2/2.d0)
        p_r = gamma1 * (ql(3,i  ) - rho_r*u_r**2/2.d0)
        
        ! Sound speed
        a_l = dsqrt(gamma*p_l/rho_l)
        a_r = dsqrt(gamma*p_r/rho_r)

        !Total energy
        E_l = qr(3,i-1)
        E_r = ql(3,i  )
     	
     	!Checking for dry states

     
 	    d_l=u_l+2*a_l/gamma1
 	    d_r=u_r-2*a_r/gamma1
        !First we check if there is a vacuum (or close to vacuum) state
        if (p_l<1.e-30 .AND. p_r<1.e-30) then !If both states are vacuum do nothing
            p_bar=0.
            h_bar=0.
            rho_bar=0.
            switch=0
        else if (p_r<tol .OR. rho_r<tol) then
            ss3=0.
            cs11=u_l-a_l
            cs12=d_l
            rho_ml=rho_l
            rho_mr=rho_r
            if (0<=cs11) then
                p_bar=p_l
                u_bar=u_l
                rho_bar=rho_l
                
                switch=5 !1-non transonic rarefaction and 3-shock 
            else if (cs12<=0.) then
                p_bar=p_r
                u_bar=u_r
                rho_bar=rho_r
                switch=5  !1-non transonic rarefaction and 3-shock
            else
                rho_bar=rho_l*(2/(gamma+1)+gamma1*u_l/((gamma+1)*a_l))**(2/gamma1)
                u_bar=2*(a_l+gamma1*u_l/2.)/(gamma+1)
                p_bar=p_l*(2/(gamma+1)+gamma1*u_l/((gamma+1)*a_l))**(2/gamma1)
                switch= 5 !1-transonic rarefaction and 3-shock
            end if
            
        else if (p_l<tol .OR. rho_l<tol) then
            ss1=0.
            cs31=d_r
            cs32=u_r+a_r
            rho_ml=rho_l
            rho_mr=rho_r
            if(0.<=cs31) then
                p_bar=p_l
                u_bar=u_l
                rho_bar=rho_l
                switch= 2   !1-shock and 3-nontransonic rarefaction
            else if (cs32<=0.) then
                p_bar=p_r
                u_bar=u_r
                rho_bar=rho_r
                switch=2 !1-shock and 3-nontransonic rarefaction
            else
                p_bar=p_r*(2/(gamma+1)-gamma1*u_r/((gamma+1)*a_r))**(2*gamma/gamma1)
                rho_bar=rho_r*(2/(gamma+1)-gamma1*u_r/((gamma+1)*a_r))**(2/gamma1)
                u_bar=2*(-a_r+gamma1*u_r/2.d0)/(gamma+1)
                switch=7 !1-shock and 3-transonic rarefaction
            end if
            
        else if (d_l<=d_r) then
            cs11=u_l-a_l
            cs12=d_l
            cs31=u_r-a_r
            cs32=d_r
            rho_ml=rho_l
            rho_mr=rho_r
            if (0<=d_l) then
                rho_bar=rho_l*(2/(gamma+1)+gamma1*u_l/((gamma+1)*a_l))**(2/gamma1)
                u_bar=2*(a_l+gamma1*u_l/2.)/(gamma+1)
                p_bar=p_l*(2/(gamma+1)+gamma1*u_l/((gamma+1)*a_l))**(2/gamma1)
            else if (d_r<=0) then
                p_bar=p_r*(2/(gamma+1)-gamma1*u_r/((gamma+1)*a_r))**(2*gamma/gamma1)
                rho_bar=rho_r*(2/(gamma+1)-gamma1*u_r/((gamma+1)*a_r))**(2/gamma1)
                u_bar=2*(-a_r+gamma1*u_r/2.d0)/(gamma+1)
            else
                p_bar=0.
                u_bar=0.
                rho_bar=0.
            end if
        else            
            ! If there are not dry states, we use the solver for wet states
            p_m=pStarNonVacuum(gamma,p_l,u_l,rho_l,p_r,u_r,rho_r,a_l,a_r,tol) 
            u_m=uStar(gamma,p_m,p_l,u_l,rho_l,p_r,u_r,rho_r,a_l,a_r)


            !Finding w_bar with non vacuum sates
            if (0<=u_m) then !If 0=x/t is at left of contact wave
                if (p_m > p_l) then   !If 1-wave is shock
                    !Left shock
                    ss1=u_l-a_l*dsqrt((gamma+1)*p_m/(2*gamma*p_l)+0.5d0*gamma1/gamma)
                    rho_ml=rho_l*(p_m/p_l+gamma1/(gamma+1))/(gamma1*p_m/((gamma+1)*p_l)+1)
                    if (0<ss1) then  !If 1-shock is moving right, wbar=wL
                        p_bar=p_l
                        u_bar=u_l
                        rho_bar=rho_l 
                    else    !If 1-shock is moving left, wbar=w*L
                        p_bar=p_m
                        u_bar=u_m
                        rho_bar=rho_ml
                    end if
                    if(p_m>p_r) then !If 3-wave is a shock
                        ss3=u_r+a_r*dsqrt((gamma+1)*p_m/(2*gamma*p_r)+0.5d0*gamma1/gamma)
                        rho_mr=rho_r*(p_m/p_r+gamma1/(gamma+1))/(gamma1*p_m/((gamma+1)*p_r)+1)
                        switch=1 !two shocks
                    else !If 3-wave is a rarefaction
                        cs31=u_m+a_r*(p_m/p_r)**(0.5d0*gamma1/gamma)
                        cs32=u_r+a_r
                        rho_mr=rho_r*(p_m/p_r)**(1/gamma)
                        switch=2 !1-shock and 3-nontransonic rarefaction
                    end if
                
                else  !If 1-wave is rarefaction
                    cs11=u_l-a_l
                    cs12=u_m-a_l*(p_m/p_l)**(0.5d0*gamma1/gamma)
                    rho_ml=rho_l*(p_m/p_l)**(1/gamma)
                    if (cs11>0) then !Wbar=WL
                        p_bar=p_l
                        u_bar=u_l
                        rho_bar=rho_l
                        if(p_m>p_r) then !If 3-wave is a shock
                            ss3=u_r+a_r*dsqrt((gamma+1)*p_m/(2*gamma*p_r)+0.5d0*gamma1/gamma)
                            rho_mr=rho_r*(p_m/p_r+gamma1/(gamma+1))/(gamma1*p_m/((gamma+1)*p_r)+1)
                            switch=3 !1-nontransonic rarefaction and 3-shock
                        else !If 3-wave is a rarefaction
                            cs31=u_m+a_r*(p_m/p_r)**(0.5d0*gamma1/gamma)
                            cs32=u_r+a_r
                            rho_mr=rho_r*(p_m/p_r)**(1/gamma)
                            switch=4 !two nontransonic rarefactions
                        end if
                    else if (cs12<0) then !Wbar=W*L
                        p_bar=p_m
                        u_bar=u_m
                        rho_bar=rho_ml
                        if(p_m>p_r) then !If 3-wave is a shock
                            ss3=u_r+a_r*dsqrt((gamma+1)*p_m/(2*gamma*p_r)+0.5d0*gamma1/gamma)
                            rho_mr=rho_r*(p_m/p_r+gamma1/(gamma+1))/(gamma1*p_m/((gamma+1)*p_r)+1)
                            switch=3 !1-nontransonic rarefaction and 3-shock
                        else !If 3-wave is a rarefaction
                            cs31=u_m+a_r*(p_m/p_r)**(0.5d0*gamma1/gamma)
                            cs32=u_r+a_r
                            rho_mr=rho_r*(p_m/p_r)**(1/gamma)
                            switch=4 !two nontransonic rarefactions
                        end if
                    else !transonic 1-rarefaction, x/t=0 inside fan
                        p_bar=p_l*(2/(gamma+1)+gamma1*u_l/((gamma+1)*a_l))**(2*gamma/gamma1)
                        rho_bar=rho_l*(2/(gamma+1)+gamma1*u_l/((gamma+1)*a_l))**(2/gamma1)
                        u_bar=2*(a_l+gamma1*u_l/2.d0)/(gamma+1)
                        if(p_m>p_r) then !If 3-wave is a shock
                            ss3=u_r+a_r*dsqrt((gamma+1)*p_m/(2*gamma*p_r)+0.5d0*gamma1/gamma)
                            rho_mr=rho_r*(p_m/p_r+gamma1/(gamma+1))/(gamma1*p_m/((gamma+1)*p_r)+1)
                            switch=5 !1-transonic rarefaction and 3-shock
                        else !If 3-wave is a rarefaction
                            cs31=u_m+a_r*(p_m/p_r)**(0.5d0*gamma1/gamma)
                            cs32=u_r+a_r
                            rho_mr=rho_r*(p_m/p_r)**(1/gamma)
                            switch=6 !1-transonic and 3-nontransonic rarefactions
                        end if
                    end if
                    
                    
                end if
                
            else !Contact wave is going to the left
                if (p_m>p_r) then!3-wave is a shock
                    ss3=u_r+a_r*dsqrt((gamma+1)*p_m/(2*gamma*p_r)+0.5d0*gamma1/gamma)
                    rho_mr=rho_r*(p_m/p_r+gamma1/(gamma+1))/(gamma1*p_m/((gamma+1)*p_r)+1)
                    if (ss3>0) then !if 3-shock is going to the right
                        p_bar=p_m
                        u_bar=u_m
                        rho_bar=rho_mr
                    else !If 3-shock is going to the left
                        p_bar=p_r
                        u_bar=u_r
                        rho_bar=rho_r
                    end if
                    if (p_m > p_l) then   !If 1-wave is shock
                    !Left shock
                        ss1=u_l-a_l*dsqrt((gamma+1)*p_m/(2*gamma*p_l)+0.5d0*gamma1/gamma)
                        rho_ml=rho_l*(p_m/p_l+gamma1/(gamma+1))/(gamma1*p_m/((gamma+1)*p_l)+1)
                        switch=1 !two shocks
                    else
                        cs11=u_l-a_l
                        cs12=u_m-a_l*(p_m/p_l)**(0.5d0*gamma1/gamma)
                        rho_ml=rho_l*(p_m/p_l)**(1/gamma)
                        switch=3 !1-nontransonic rarefaction and 3-shock 
                    end if 
                else !if 3-wave is a rarefaction
                    cs31=u_m+a_r*(p_m/p_r)**(0.5d0*gamma1/gamma)
                    cs32=u_r+a_r
                    rho_mr=rho_r*(p_m/p_r)**(1/gamma)
                    if (cs31>0) then !3-rarefaction going to the right
                        p_bar=p_m
                        u_bar=u_m
                        rho_bar=rho_mr
                        if (p_m > p_l) then   !If 1-wave is shock
                            !Left shock
                            ss1=u_l-a_l*dsqrt((gamma+1)*p_m/(2*gamma*p_l)+0.5d0*gamma1/gamma)
                            rho_ml=rho_l*(p_m/p_l+gamma1/(gamma+1))/(gamma1*p_m/((gamma+1)*p_l)+1)
                            switch=2 !1-shock and 3-nontransonic rarefaction
                        else
                            cs11=u_l-a_l
                            cs12=u_m-a_l*(p_m/p_l)**(0.5d0*gamma1/gamma)
                            rho_ml=rho_l*(p_m/p_l)**(1/gamma)
                            switch=4 !two nontransonic rarefactions
                        end if 
                    else if (cs32<0) then !3-rarefaction going to the left
                        p_bar=p_r
                        u_bar=u_r
                        rho_bar=rho_r
                        if (p_m > p_l) then   !If 1-wave is shock
                            !Left shock
                            ss1=u_l-a_l*dsqrt((gamma+1)*p_m/(2*gamma*p_l)+0.5d0*gamma1/gamma)
                            rho_ml=rho_l*(p_m/p_l+gamma1/(gamma+1))/(gamma1*p_m/((gamma+1)*p_l)+1)
                            switch=2 !1-shock and 3-nontransonic rarefaction
                        else
                            cs11=u_l-a_l
                            cs12=u_m-a_l*(p_m/p_l)**(0.5d0*gamma1/gamma)
                            rho_ml=rho_l*(p_m/p_l)**(1/gamma)
                            switch=4 !two nontransonic rarefactions
                        end if 
                    else   !transonic 3-rarefaction
                        p_bar=p_r*(2/(gamma+1)-gamma1*u_r/((gamma+1)*a_r))**(2*gamma/gamma1)
                        rho_bar=rho_r*(2/(gamma+1)-gamma1*u_r/((gamma+1)*a_r))**(2/gamma1)
                        u_bar=2*(-a_r+gamma1*u_r/2.d0)/(gamma+1)
                        if (p_m > p_l) then   !If 1-wave is shock
                            !Left shock
                            ss1=u_l-a_l*dsqrt((gamma+1)*p_m/(2*gamma*p_l)+0.5d0*gamma1/gamma)
                            rho_ml=rho_l*(p_m/p_l+gamma1/(gamma+1))/(gamma1*p_m/((gamma+1)*p_l)+1)
                            switch=7 !1-shock and 3-transonic rarefaction
                        else
                            cs11=u_l-a_l
                            cs12=u_m-a_l*(p_m/p_l)**(0.5d0*gamma1/gamma)
                            rho_ml=rho_l*(p_m/p_l)**(1/gamma)
                            switch=8 !1-nontransonic and 3-transonic rarefactions
                        end if 
                    end if
                end if
                 
            end if 
        end if                
                        
                        
        !Obtaining q_bar if not done yet, computing speeds and waves
        !If switch < 5 there are not transonic rarefactions, q bar has not been found
        !If switch >= 5 there is a transonic rarefaction, q bar has been found
                
        select case (switch)
            case (0) !Do nothing
                s(1,i)=0.
                s(2,i)=0.
                s(3,i)=0.
                
                wave(1,1,i)=0.
                wave(2,1,i)=0.
                wave(3,1,i)=0.
                
                wave(1,2,i)=0.
                wave(2,2,i)=0.
                wave(3,2,i)=0.
                
                wave(1,3,i)=0.
                wave(2,3,i)=0.
                wave(3,3,i)=0.
            case (1) !both waves are shocks
                
                s(1,i)=ss1
                s(2,i)=u_m 
                s(3,i)=ss3
                
                wave(1,1,i)=rho_ml-qr(1,i-1)
                wave(2,1,i)=rho_ml*u_m-qr(2,i-1)
                wave(3,1,i)=rho_ml*0.5d0*u_m**2+p_m/gamma1-qr(3,i-1)
                
                wave(1,2,i)=rho_mr-rho_ml
                wave(2,2,i)=rho_mr*u_m-rho_ml*u_m
                wave(3,2,i)=rho_mr*0.5d0*u_m**2-rho_ml*0.5d0*u_m**2
                
                wave(1,3,i)=ql(1,i  )-rho_mr
                wave(2,3,i)=ql(2,i  )-rho_mr*u_m
                wave(3,3,i)=ql(3,i  )-rho_mr*0.5d0*u_m**2-p_m/gamma1
                
            case (2)  !1-shock and 3-nontransonic rarefaction
            
                s(1,i)=ss1
                s(2,i)=u_m 
                s(3,i)=cs32
                
                wave(1,1,i)=rho_ml-qr(1,i-1)
                wave(2,1,i)=rho_ml*u_m-qr(2,i-1)
                wave(3,1,i)=rho_ml*0.5d0*u_m**2+p_m/gamma1-qr(3,i-1)
                
                wave(1,2,i)=rho_mr-rho_ml
                wave(2,2,i)=rho_mr*u_m-rho_ml*u_m
                wave(3,2,i)=rho_mr*0.5d0*u_m**2-rho_ml*0.5d0*u_m**2
                
                wave(1,3,i)=(rho_r*u_r-rho_mr*u_m)/s(3,i)
                wave(2,3,i)=(rho_r*u_r**2+p_r-rho_mr*u_m**2-p_m)/s(3,i)
                wave(3,1,i)=(u_r*(E_r+p_r)-u_m*(0.5d0*rho_mr*u_m**2+p_m/gamma1+p_m))/s(3,i)
                
                
            case (3)   !1-non transonic rarefaction, 3-shock
                s(1,i)=cs11
                s(2,i)=u_m 
                s(3,i)=ss3
                
                wave(1,1,i)=(rho_ml*u_m-rho_l*u_l)/s(1,i)
                wave(2,1,i)=(rho_ml*u_m**2+p_m-rho_l*u_l**2-p_l)/s(1,i)
                wave(3,1,i)=(u_m*(0.5d0*rho_ml*u_m**2+p_m/gamma1+p_m)-&
                            u_l*(E_l+p_l))/s(1,i)
                            
                wave(1,2,i)=rho_mr-rho_ml
                wave(2,2,i)=rho_mr*u_m-rho_ml*u_m
                wave(3,2,i)=rho_mr*0.5d0*u_m**2-rho_ml*0.5d0*u_m**2
                
                wave(1,3,i)=ql(1,i  )-rho_mr
                wave(2,3,i)=ql(2,i  )-rho_mr*u_m
                wave(3,3,i)=ql(3,i  )-rho_mr*0.5d0*u_m**2-p_m/gamma1
                
            case (4)   !two nontransonic rarefactions
                
                s(1,i)=cs11
                s(2,i)=u_m 
                s(3,i)=cs32
                
                wave(1,1,i)=(rho_ml*u_m-rho_l*u_l)/s(1,i)
                wave(2,1,i)=(rho_ml*u_m**2+p_m-rho_l*u_l**2-p_l)/s(1,i)
                wave(3,1,i)=(u_m*(0.5d0*rho_ml*u_m**2+p_m/gamma1+p_m)-&
                            u_l*(E_l+p_l))/s(1,i)
                
                wave(1,2,i)=rho_mr-rho_ml
                wave(2,2,i)=rho_mr*u_m-rho_ml*u_m
                wave(3,2,i)=rho_mr*0.5d0*u_m**2-rho_ml*0.5d0*u_m**2
                            
                wave(1,3,i)=(rho_r*u_r-rho_mr*u_m)/s(3,i)
                wave(2,3,i)=(rho_r*u_r**2+p_r-rho_mr*u_m**2-p_m)/s(3,i)
                wave(3,1,i)=(u_r*(E_r+p_r)-u_m*(0.5d0*rho_mr*u_m**2+p_m/gamma1+p_m))/s(3,i)
                
        end select
        
        !Computing fluctuations now, since we will use them to compute waves for
        !transonic rarefactions

        amdq(1,i)=rho_bar*u_bar-rho_l*u_l
        amdq(2,i)=rho_bar*u_bar**2+p_bar-rho_l*u_l**2-p_l
        amdq(3,i)=u_bar*(0.5d0*rho_bar*u_bar**2+p_bar/gamma1+p_bar)-u_l*(E_l+p_l)
                            
                            
        apdq(1,i)=rho_r*u_r-rho_bar*u_bar
        apdq(2,i)=rho_r*u_r**2+p_r-rho_bar*u_bar**2-p_bar
        apdq(3,i)=u_r*(E_r+p_r)-u_bar*(0.5d0*rho_bar*u_bar**2+p_bar/gamma1+p_bar)
        
        !We continue computing waves and speeds
        select case (switch)
            case(5) !1-transonic, 3-shock
                s(1,i)=cs11
                s(2,i)=0.d0
                s(3,i)=ss3
                
                wave(1,1,i)=amdq(1,i)/s(1,i)
                wave(2,1,i)=amdq(2,i)/s(1,i)
                wave(3,1,i)=amdq(3,i)/s(1,i)
                
                wave(1,2,i)=0.d0
                wave(2,2,i)=0.d0
                wave(3,2,i)=0.d0
                
                !wave(1,3,i)=ql(1,i  )-rho_bar
                !wave(2,3,i)=ql(2,i  )-rho_bar*u_bar
                !wave(3,3,i)=ql(3,i  )-0.5d0*rho_bar*u_bar**2-p_bar/gamma1
                
                wave(1,3,i)=apdq(1,i)/s(3,i)
                wave(2,3,i)=apdq(2,i)/s(3,i)
                wave(3,3,i)=apdq(3,i)/s(3,i)
                
            case(6) !1-transonic, 3-nontransonic rarefactions
                s(1,i)=cs11
                s(2,i)=0.d0
                s(3,i)=cs32
                
                wave(1,1,i)=amdq(1,i)/s(1,i)
                wave(2,1,i)=amdq(2,i)/s(1,i)
                wave(3,1,i)=amdq(3,i)/s(1,i)
                
                wave(1,2,i)=0.d0
                wave(2,2,i)=0.d0
                wave(3,2,i)=0.d0
                
                wave(1,3,i)=apdq(1,i)/s(3,i)
                wave(2,3,i)=apdq(2,i)/s(3,i)
                wave(3,3,i)=apdq(3,i)/s(3,i)
            
            case(7) !1-shock 3-transonic raref
                s(1,i)=ss1
                s(2,i)=0.d0
                s(3,i)=cs32
                
                !wave(1,1,i)=rho_bar-qr(1,i-1)
                !wave(2,1,i)=rho_bar*u_bar-qr(2,i-1)
                !wave(3,1,i)=0.5d0*rho_bar*u_bar**2+p_bar/gamma1-qr(3,i-1)
                
                wave(1,1,i)=amdq(1,i)/s(1,i)
                wave(2,1,i)=amdq(2,i)/s(1,i)
                wave(3,1,i)=amdq(3,i)/s(1,i)
                
                wave(1,2,i)=0.d0
                wave(2,2,i)=0.d0
                wave(3,2,i)=0.d0
                
                wave(1,3,i)=apdq(1,i)/s(3,i)
                wave(2,3,i)=apdq(2,i)/s(3,i)
                wave(3,3,i)=apdq(3,i)/s(3,i)
            case(8) !1-nontransonic, 3-transonic
                s(1,i)=cs11
                s(2,i)=0.d0
                s(3,i)=cs32
                
                wave(1,1,i)=amdq(1,i)/s(1,i)
                wave(2,1,i)=amdq(2,i)/s(1,i)
                wave(3,1,i)=amdq(3,i)/s(1,i)
                
                wave(1,2,i)=0.d0
                wave(2,2,i)=0.d0
                wave(3,2,i)=0.d0
                
                wave(1,3,i)=apdq(1,i)/s(3,i)
                wave(2,3,i)=apdq(2,i)/s(3,i)
                wave(3,3,i)=apdq(3,i)/s(3,i)
        end select
!     
    end do 

end subroutine rp1

double precision function TwoRarefactionInitialGuess(gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol) 
  !This the exact solution based on having two rarefactions
  implicit double precision (a-h,o-z)
  gamma1=gamma-1.d0
  TwoRarefactionInitialGuess = ((aL+aR-0.5d0*gamma1*(uR-uL))/(aL/pL**(gamma1/(2*gamma))  &
                               +aR/pR**(gamma1/(2*gamma))))**(2*gamma/gamma1)

END function TwoRarefactionInitialGuess



double precision function pStarNonVacuum(gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol) 

    implicit none

    double precision, intent(in) :: gamma,pL,uL,rhoL,pR,uR,rhoR,aL,aR,tol
    integer :: iterations

    double precision :: CL,CR,BL,BR
    double precision ::pMin, phiMin, TwoRarefactionInitialGuess, phi, phip, extra, ustar

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
        ! We have one rarefaction and one shock or two shocks
        ppv=max(pMin,0.5d0*(pL+pR)-0.125d0*(uR-uL)*(rhoL+rhoR)*(aL+aR))

        p = (gs(gamma, ppv, CL, BL)*pL+gs(gamma, ppv, CR, BR)*pR-(uR-uL))/&
            (gs(gamma,ppv,CL,BL)+gs(gamma,ppv,CR,BR)) 
    END if

    if (iterate==1) then 
        p = max(pMin,p-phi(gamma,p,pL,uL,rhoL,pR,uR,rhoR,aL,aR,CL,CR,BL,BR)/&
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
END function pStarNonVacuum

double precision function gs(gamma, p, CS, BS)
   implicit double precision (a-h,o-z)
   gs=dsqrt(CS/(p+BS))
end function gs

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
