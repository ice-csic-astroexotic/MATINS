c       
c       This routine takes the velocity field (FB / nc) and 
c       Calculates 
c       1. the chemical equilibrium deviation 
c       Delta mu = mu_p + mu_e - mu_n
c       2. The ambipolar diffusion velocity 
c       Vamb = coef_amb * [ V - grad(Delta mu) ] 
c
c       Note: It needs the following coefficients in a common block
c       coef_a(m) =   1d0 / aa_tmp(m)**2
c       coef_b(m) =   1d0 / bb_tmp(m)         
c       coef_amb(m) = omega_p * tau_tmp(m)

c       Note: this routine uses the subroutine 
c       "div" to determine the divergence of V (Source)
c       "gradient" 
c       to determine the grad (Delta mu)

c       -------------------------------
c	subroutine operator_spherical
      subroutine ambipolar(lcore,Vin,Vamb)
c       -------------------------------
      use structure, only: rho,xh,ye,yn,yp,aa,zz
      implicit none

!-----------------------------------------------------------------------
! Dimensions of the grid.
!-----------------------------------------------------------------------
      include '../decl/dim2.h'
      
      integer nx,nz, lcore
      real*8 rd(lmax),td(kmax-1)

      integer m,n, i,j
      real*8 dr, dt, coef1, dfdr, dfdt
      real*8 rhomed, ynmed, ypmed

      real*8, dimension(3,lmax,kmax-1) :: grad_Phi,Vin,Vamb
      real*8, dimension(lmax) :: lambda,lambda_DU,tau_pn,beta,alpha
      real*8, dimension(lmax) :: coef_a,coef_b,coef_amb
C
C     Arrays to be passed to TRDIG are dimensioned (lcore x kmax)
C
      real*8, dimension(lcore,kmax-1) :: Phi,Src,div_V
      real*8 a(lcore,kmax-1,5)
      real*8 rho_nuc, mp, nbaryon, omega_p, t8

c       a1 * f(m,n-1) + a2 * f(m-1,n) + a3 * f(m,n) + a4 * f(m+1,n) + a5 * f(m,n+1) = Src

c                               f(m-1,n-1)
c                               f(m  ,n-1)   
c	                        f(m+1,n-1) 
c	..   a1    ..           f(m-1,n)
c	a2   a3    a4           f(m  ,n)
c       ..   a5    ..           f(m+1,n)
c---------------------------------------------------------------
c       1. Solve the Modified Poisson equation 

c	opL (Phi) = Src,  where
c       opL:  nabla^2 + coef_b * d_dr - coef_a 
c
c       and coef_a = 1 / a**2
c
c       and coef_b = 1 / b
c       1/b = dln(beta)_dr
c
c       The solution provide the chemical potential 
c       deviation:  Phi = Delta mu
c-----------------------------------------------------------
      nx=lcore
      nz=kmax-1

      do m=1,lmax
        rd(m)=rb(2*m-1)
      enddo
      do n=1,nz
        td(n)=theta(2*n)
      enddo

      dr = rd(2)-rd(1)
      dt = td(2)-td(1)

      t8=10.d0
      rho_nuc=2.8d14
      mp = 1.672d-24  ! proton mass
      omega_p = 9.6043d15   ! Proton Larmor frequency

      CALL ambipolar_sf(t8,lambda,lambda_DU,tau_pn)

      do m=1,nx
        rhomed=rho(2*m)
	      ynmed=yn(2*m)
	      ypmed=yp(2*m)
c        lambda(m) = 5d27*t8**6*(rhomed(m)/rho_nuc)**(2d0/3d0)
c        lambda_DU(m)=3.5d36*t8**4*(rhomed(m)/rho_nuc)**(1d0/3d0)
c        tau_pn(m)=(rhomed(m)/rho_nuc)**(1d0/3d0)/(4.7d16*t8**2)
        beta(m)=ynmed*rhomed*ypmed*tau_pn(m)/mp**2
c        beta(m)=1.02d34*(rhomed(m)/rho_nuc)**(4d0/3d0)/t8**2
      enddo

C-----------------------------------------------------------
C-----------------------------------------------------------

      m=1
      coef_b(m) =  (beta(m+1)-beta(m))/(beta(m)*(rd(m+1)-rd(m)))
      do m=2,nx-1
        coef_b(m) = (beta(m+1)-beta(m-1))/(beta(m)*(rd(m+1)-rd(m-1)))
      enddo
      m=nx
      coef_b(m) =  (beta(m)-beta(m-1))/(beta(m)*(rd(m)-rd(m-1)))

      do m=1,nx
        ynmed=yn(2*m)
        coef_a(m) = lambda(m)/(ynmed*beta(m))*1.d10  !in km**-2
        coef_amb(m) = omega_p*tau_pn(m)
      enddo

c       ---------------------------------
c       Matrix in the internal points
c       --------------------------------

      a=0d0

      do n=1,nz		! theta  loop
      do m=1,nx		! radial loop

        coef1=0.5d0*dcos(td(n))/(rd(m)**2*dsin(td(n))*dt)
        alpha(m)  = 1d0 + 0.25d0 * (dr/rd(m)) **2

        a(m,n,1) = 1d0 / ( rd(m) * dt )**2 - coef1

        a(m,n,2) = alpha(m)/dr**2 - 1d0/(rd(m)*dr)
     &                    - 0.5d0*coef_b(m)/dr

        a(m,n,3) = - 2d0*alpha(m)/dr**2 
     &                    - 2d0/(rd(m)*dt)**2
     &                    - coef_a(m) 

        a(m,n,4) = alpha(m)/dr**2 + 1d0/(rd(m)*dr)
     &                    + 0.5d0*coef_b(m)/dr

        a(m,n,5) = 1d0/(rd(m)*dt)**2 + coef1

      enddo
      enddo

c       --------------------------------------------------
c       Determine the source of the Poisson-like equation
c       opL (Phi) = Srs = div(Vin)
c       -------------------------------------------------
      div_V = 0.d0

      do n=1,nz               ! theta  loop
      do m=1,nx            ! radial loop
       if (m.eq.1) then
         dfdr = -0.5*(3.0*Vin(1,m,n)
     &   -4.0*Vin(1,m+1,n)+Vin(1,m+2,n) )/dr 
       elseif (m.lt.nx) then
         dfdr = 0.5*(Vin(1,m+1,n)-Vin(1,m-1,n))/dr
       else
         dfdr = 0.5*(3.0*Vin(1,m,n)
     &       -4.0*Vin(1,m-1,n)+Vin(1,m-2,n) )/dr
       endif
       if (n.eq.1) then
         dfdt = -0.5*(3.0*Vin(2,m,n)
     &      -4.0*Vin(2,m,n+1)+Vin(2,m,n+2) )/dt
       elseif (n.lt.nz) then
         dfdt = 0.5*(Vin(2,m,n+1)-Vin(2,m,n-1))/dt
       else
         dfdt = 0.5*(3.0*Vin(2,m,n)
     &      -4.0*Vin(2,m,n-1)+Vin(2,m,n-2))/dt
       endif

       div_V(m,n) = dfdr  + 2d0*Vin(1,m,n)/rd(m)
     &    + dfdt/rd(m)
     &    + dcos(td(n))/dsin(td(n))/rd(m)*Vin(2,m,n)
      enddo
      enddo

c	Src = div(V) + 1/b * Vr
      do m=1,nx
      do n=1,nz
        Src(m,n) = div_V(m,n) + coef_b(m)*Vin(1,m,n)
      enddo
      enddo

c       ---------------------------------
c       Matrix at the boundaries
c       --------------------------------

c       1) BC at theta = Theta_in (n=1)
c       Neumann (df_dtheta)
	   n=1
	   do m=2,nx-1
	      a(m,n,3) = a(m,n,3) + a(m,n,1) 
	      a(m,n,1) = 0d0
	   enddo

c       2) BC at theta = Theta_out (n=nz)
c       Neumann (df_dtheta)
	   n=nz
	   do m=2,nx-1
	      a(m,n,3) = a(m,n,3) + a(m,n,5) 
	      a(m,n,5) = 0d0
	   enddo

c       3) BC at r=rin (m=1).
c       Neumann (df_dr)
	   m=1
	   do n=2,nz-1
	      a(m,n,4) = a(m,n,4) + a(m,n,2) 
	      a(m,n,2) = 0d0
	   enddo

c       4) BC at r=out (m=nx).
c       Neumann (df_dr)
	   m=nx
	   do n=2,nz-1
              Src(m,n) = Src(m,n) - 2.d0*Vin(1,m,n)*dr*a(m,n,4)
	      a(m,n,2) = a(m,n,2) + a(m,n,4) 
	      a(m,n,4) = 0d0
	   enddo

	   m=1
	   n=1
	   a(m,n,3) = a(m,n,3) + a(m,n,1) + a(m,n,2)
	   a(m,n,1) = 0d0
	   a(m,n,2) = 0d0

	   m=1
	   n=nz
	   a(m,n,3) = a(m,n,3) + a(m,n,2) + a(m,n,5)
	   a(m,n,2) = 0d0
	   a(m,n,5) = 0d0

	   m=nx
	   n=1
	   a(m,n,3) = a(m,n,3) + a(m,n,1) + a(m,n,4)
	   a(m,n,1) = 0d0
	   a(m,n,4) = 0d0

	   m=nx
	   n=nz
	   a(m,n,3) = a(m,n,3) + a(m,n,4) + a(m,n,5)
	   a(m,n,4) = 0d0
	   a(m,n,5) = 0d0

!-----------------------------------------------------------------------
! Solution of the linear system.
        Phi=0.d0
        call trdig(nx,nz,Phi,a,Src)

        Phi=Phi-Phi(1,1)

c-----------------------------------------------------------------------
c       Calculation of the Ambipolar velocity
c-----------------------------------------------------------------------

      grad_Phi = 0.d0
c       1. Gradient of Phi where (Phi = Delta mu) 
      do n=1,nz               ! theta  loop
      do m=1,nx            ! radial loop
       if (m.eq.1) then
         dfdr = -0.5*(3.0*Phi(m,n)-4.0*Phi(m+1,n)+Phi(m+2,n) )/dr
       elseif (m.lt.nx) then
         dfdr = 0.5*(Phi(m+1,n)-Phi(m-1,n))/dr
       else
         dfdr = 0.5*(3.0*Phi(m,n)-4.0*Phi(m-1,n)+Phi(m-2,n) )/dr
       endif
       if (n.eq.1) then
         dfdt = -0.5*(3.0*Phi(m,n)-4.0*Phi(m,n+1)+Phi(m,n+2) )/dt
       elseif (n.lt.nz) then
         dfdt = 0.5*(Phi(m,n+1)-Phi(m,n-1))/dt
       else
         dfdt = 0.5*(3.0*Phi(m,n)-4.0*Phi(m,n-1)+Phi(m,n-2) )/dt
       endif

       grad_Phi(1,m,n) = dfdr 
       grad_Phi(2,m,n) = dfdt/rd(m)
      enddo
      enddo

c       2. Ambipolar velocity Vamb = coef_amb * [ V - grad(Delta mu) ] 
      Vamb = 0.d0
      Vamb = Vin !- grad_Phi

      do m=1,nx
        Vamb(1,m,:)=Vamb(1,m,:)*coef_amb(m)
        Vamb(2,m,:)=Vamb(2,m,:)*coef_amb(m)
        Vamb(3,m,:)=Vamb(3,m,:)*coef_amb(m)
      enddo

!---------------------------------------
      return
      end


