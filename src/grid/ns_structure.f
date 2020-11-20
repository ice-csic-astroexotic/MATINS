!!
!!        This programs solves the ODEs (TOV equations) to obtain 
!!       the structure of a neutron star, using a tabulated EOS
!! TODO: - Provide radial grid suitable for different central pressure (different masses)
!!
      subroutine solve_TOV_structure(nradial, radius, xmass, pres, g14)

      use input_params, only: use_relativistic_grid
      use grid, only: benu, belam
      use constants, only: PI, c2dg, c4dg
      use legpol, only: get_rel_correction
      implicit none 

      integer nradial
      real*8 radius(nradial), pres(nradial), xmass(nradial)

      integer i,nok,nbad
      external derivs, derivs2
      integer, parameter :: nvar=3
      real*8 y(nvar),y0(nvar), x1,x2,h1,hh
      real*8 rhoeos, rho0, z, a, xn, xh, pmin
      real*8 g14, double_compactness
      real*8 pc, rhoc, pcgs, nu0

      pcgs = 1.37d35
      write(*,*) '[STRUCTURE]: Central pressure = ', pcgs
      pmin=1.d12
      
      call geteost(pcgs,rhoeos,rho0,z,a,xn,xh) 
      PC = pcgs/c4dg  ! cm**-2
      rhoc = rhoeos/c2dg  ! cm**-2

!----------------------------------------------------
!  Integrate once down to very low pressure 
!   to obtain radius, mass, nu0, etc.
!----------------------------------------------------
      hh = 1.d-2*PC
      x1 = DLOG(PC-hh)
      x2 = DLOG(pmin/c4dg)
!----------------------------------------------------              
!      Choose where to stop (crust)
!----------------------------------------------------              
!      x2 = DLOG(4.1d-21)   ! Pre in cm^-2 para rho=3e10 g/cc
!      x2 = DLOG(1.7d-21)   ! Pre in cm^-2 para rho=1e10 g/cc
!      x2 = DLOG(4.d-23)   ! Pre in cm^-2 para rho=1e9 g/cc
!      x2 = DLOG(3.d-24)   ! Pre in cm^-2 para rho=1e8 g/cc
!      x2 = DLOG(1.3d-26)   ! Pre in cm^-2 para rho=5e6 g/cc
!      x2 = DLOG(4.d-27)   ! Pre in cm^-2 para rho=1e6 g/cc
!      x2 = DLOG(3.d-29)   ! Pre in cm^-2 para rho=5e6 g/cc

      Y0(1) = DSQRT((2.d0*hh)/((rhoc+PC)*(rhoc/3.d0+PC)*4.d0*pi))
      Y0(2) = 4.d0*pi*rhoc*Y0(1)**3/3.d0
      Y0(3) = 0.d0

      Y = Y0

      h1=1.d-2*(x2-x1)

      call odeint(y,nvar,x1,x2,h1,nok,nbad,derivs)
      
      nu0 = -y(3)+dlog(1.d0-2.d0*y(2)/y(1))
 
!----------------------------------------------------
!       Second integration, by steps
!----------------------------------------------------
      x1 = 1.d-2*radius(1)*1.d5
      Y0(1) = PC-((rhoc+PC)*(rhoc/3.d0+PC)*2.d0*pi)*x1**2
      Y0(2) = 4.d0*pi*rhoc*x1**3/3.d0
      Y0(3) = nu0
      Y = Y0

      do i =1,nradial
        x2=radius(i)*1.d5
        h1=1.d-2*(x2-x1)
        call odeint(y,nvar,x1,x2,h1,nok,nbad,derivs2)
        x1=x2
        pcgs = Y(1)*c4dg
        call geteost(pcgs,rhoeos,rho0,z,a,xn,xh) 
        pres(i) = pcgs
        benu(i) = dexp(Y(3)/2.d0)
        belam(i) = 1.d0/dsqrt(1.d0-2.d0*Y(2)/x2)
        xmass(i)= Y(2)/1.4766d5     ! Solar masses
      enddo

      ! Set to one the Relativistic factors e^lam and e^nu (newtonian limit).
      if (use_relativistic_grid .eqv. .true.) then
        write(*,*) "<info>", "[STRUCTURE]: ", "Relativistic grid."
        ! Nothing specifically done, handled by default.
      else
        write(*,"(a)") "<info>", "[STRUCTURE]: ", "Newtonian grid."
        benu = 1d0
        belam = 1d0
      end if

      benu(0) = benu(1)
      belam(0) = 1.d0
      g14 = 1.d-14*Y(2)/x2*(3.d10**2/x2)/benu(nradial)  ! gravity at the surface in units of 10**14 cm/s**2
      double_compactness = 1d0 - 1d0/belam(nradial)**2 ! 2*G*M/c**2*R
      
      call get_rel_correction(double_compactness)
  

      end subroutine solve_TOV_structure


      subroutine derivs(x,y,dydx)
      use constants, only: PI, c2dg, c4dg
      implicit none
      integer nmax
      parameter (nmax=10)
      real*8 x,y(nmax),dydx(nmax)
      real*8 p,rho,z,a,xn,xh,dpdr,r,mg
      real*8 pcgs, rhoeos, rho0

      r=Y(1)
      mg=Y(2)
      p=dexp(x)
      
      pcgs = p*c4dg
      call geteost(pcgs,rhoeos,rho0,z,a,xn,xh)
      rho = rhoeos/c2dg

      dpdr = -(rho+p)*(mg+4.d0*pi*r**3*p)/(r*r-2.d0*mg*r)
      dydx(1)= p/dpdr
      dydx(2)= 4.d0*pi*r*r*rho*(p/dpdr)
      dydx(3) = -2.d0*p/(rho+p) 

      return
      end subroutine derivs

      subroutine derivs2(x,y,dydx)
      use constants, only: PI, c2dg, c4dg
      implicit none
      integer nmax
      parameter (nmax=10)
      real*8 x,y(nmax),dydx(nmax)
      real*8 p,rho,z,a,xn,xh,dpdr,r,mg
      real*8 pcgs, rhoeos, rho0

      r=x
      p=Y(1)
      mg=Y(2)
      
      pcgs = p*c4dg
      call geteost(pcgs,rhoeos,rho0,z,a,xn,xh)
      rho = rhoeos/c2dg

      dpdr = -(rho+p)*(mg+4.d0*pi*r**3*p)/(r*r-2.d0*mg*r)
      dydx(1)= dpdr
      dydx(2)= 4.d0*pi*r*r*rho
      dydx(3) = -2.d0*dpdr/(rho+p) 

      return
      end subroutine derivs2


!-----------------------------------
!    EOS (interpolation from table)
!-----------------------------------
      subroutine geteost(p,rho,n_b,z,a,xn,xh)
      implicit none
      real*8 alfa
      real*8 pt(2000),rhot(2000),nbt(2000), gam1(2000)
      real*8 zt(2000),at(2000),xnt(2000), p_temp(2000)
      real*8 p,rho, n_b, z, a, xn, xh
      integer i, j, neos

      save pt,rhot,nbt,zt,at,xnt,neos, gam1
!-------------------------------------------------
!     locate the input value of pressure in the array
!-------------------------------------------------
      !! Find the nearest pressure to "p" in the table 
      !! Uses MINLOC intrinsic function that returns the location of the
      !! minimum value in an array
      p_temp = dabs(pt - p)
      j = MINLOC(p_temp,neos)

      if( (p<pt(1)) .or. (p>pt(neos)) ) then
        write(*,'(a,i10,3es12.4)')'p is out of the table',j,p,pt(1),pt(neos)
        stop
      end if

!-------------------------------------------------
!     do logarithmic interpolation for densities
!     take the value of the nearest cell for composition 
!-------------------------------------------------
      alfa = 10.d0**(log10(p/pt(j))/gam1(j))
      rho = rhot(j)*alfa
      n_b = nbt(j)*alfa

      if (at(j)>0.d0) then
        alfa = (rho-rhot(j))/(rhot(j+1)-rhot(j-1))
        xn = xnt(j) + alfa*(xnt(j+1)-xnt(j-1))
        xn = dmax1(xn,0.d0)
        xh = 1.d0-xn
        z = zt(j) 
        a = at(j) 
      else
        xn=xnt(j) 
        xh = 0.d0
        z = 1.d0-xn
        a = 1.d0
      endif

      return

      ENTRY init_eos_tab()

      open(10,file='in/EOS_DH.tab')
      do i=1,2000
        read(10,*,end=20) nbt(i), rhot(i), pt(i), zt(i), at(i), xnt(i), gam1(i)
      enddo
 20   neos=i-1
      close(10)
      print*,'EOS TABLE LOADED: npoints=', neos

      return

      end subroutine geteost
