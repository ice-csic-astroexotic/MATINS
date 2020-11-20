!!-------------------------------------------------------------------------------
!! Module: Microphysics
!!
!!> @author
!!> Jose Pons Botella
!!> 
!!> @brief This module is responsible for all the microphysics input.
!!
!!-------------------------------------------------------------------------------
module microphysics

  use grid, only: kmax, lmax, nang, np
  use grid, only: tem0, etab, dfc2, dfc3, anis2, anis3
  use grid, only: br, bth, bphi, bm

  contains

!!-----------------------------------------------------------------------
!!   Interface with the microphysics routines
!!   Calculates diffusion coefficients, anisotropy factors and resistivity.
!!
!! Contents:
!! compute_transport
!!-----------------------------------------------------------------------
!! Cell distribution of variables.
!!
!! First true cell center marked by k=2, l=2
!! Axis in theta(i=1,nang), z(k=1,kmax)
!! (corresponding to pos.1 of k=2 and pos.1 of kmax+1)
!! (corresponding to pos.3 of k=1 and pos.3 of kmax)
!! Center of star in j=0 (not evolving variables), Rns at j=np
!! One ghost cell in angular directions (0,nang+1)
!! Two ghost cells above surface (np+1,np+2)
!!
!! Cell denoted by (k,l):
!!
!! |--------2--------|		rb(j=2l)
!! 1        0        3		rb(j=2l-1)
!! |--------4--------|		rb(j=2l-2)
!!
!! th(2k-3)  th(2k-2)  th(2k-1)
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
      subroutine compute_transport

      use constants, only: PI, CLIGHT, T_YEAR
      use grid, only: rho,xh,ye,aa,zz,tcn,tcp, yn,tccru
      use grid, only: tem0

      implicit none

! Internally used variables.
      integer i,j,k,l
      real*8 f2,kt2,kh2,kp2
      real*8 f3,kt3,kh3,kp3
      real*8 tem,taun,taup
      real*8 tcond, tcondt,tcondh,rsigma,qj,qjt,qjh
      real*8 Zimp,taucru
      real*8 CRKAPi,ksfph
!-----------------------------------------------------------------------
! Conductivity and Anisotropy factors at the vertical (2) interfaces
!-----------------------------------------------------------------------
      do l=1,lmax
      do k=1,kmax
        j = 2*l
        i = 2*k-2
        if (l == lmax) then
          tem=tem0(k,l)
        else
          tem=0.5d0*(tem0(k,l)+tem0(k,l+1))
        endif

!-----------------------------------------------------------------
        taucru=tem/dmax1(tccru(j),0.1d0*tem) 
        taun=tem/dmax1(tcn(j),0.1d0*tem)
        taup=tem/dmax1(tcp(j),0.1d0*tem)
! -----------------------------------------------------------
        IF (xh(j).eq.0) then
! -----------------------------------------------------------
!	CORE: en_cond,  No anisotropy in the core
! -----------------------------------------------------------
          call en_cond(bm(i,j),tem,rho(j),ye(j),taun,taup,tcond,etab(i,j))
          f2=0.d0
          kt2=tcond*1.d13*1.d-40 ! cgs to [10^40 erg/(s*10^8 K*km)]
          kp2=kt2
          kh2=0.d0
        ELSE
! -----------------------------------------------------------
! 	impurity parameter (input for Potekhin's routine)
! -----------------------------------------------------------
        call get_Zimp (rho(j),Zimp)
        
! -----------------------------------------------------------
!	      CRUST: potekhin2019 (condconv, condegin)         
!    electron + phonons conductivity (not SF phonons)
! -----------------------------------------------------------
       call potekhin2019(tem,rho(j),bm(i,j),aa(j),zz(j),xh(j),Zimp,ye(j), &
     & yn(j),taucru,tcond,tcondt,tcondh,rsigma,qj,qjt,qjh,CRKAPi,Ksfph)

      
        if (bm(i,j).gt.0.d0) then
          f2=(tcond/tcondt-1.d0)/bm(i,j)**2
          kh2=(tcondh/tcondt)/bm(i,j)
        else
          f2=0.d0
          kh2=0.d0
        endif
        kt2=tcondt*1.d8*1.d5*1.d-40  ! cgs to [10^40 erg/(s*10^8 K*km)]
        etab(i,j)=CLIGHT**2/(4.d0*PI*rsigma)*T_YEAR*1.D-4   ! in Km^2/Myr
        ENDIF
!-----------------------------------------------------------------------
!   Conductivity and Anisotropy factors at the lateral (3) edges.
!-----------------------------------------------------------------------
        i = 2*k-1
        j = 2*l-1
        if (k == kmax) then
          tem=tem0(k,l)
        else
          tem=0.5d0*(tem0(k,l)+tem0(k+1,l))
        endif
        taun=tem/dmax1(tcn(j),0.1d0*tem)
        taup=tem/dmax1(tcp(j),0.1d0*tem)
! -----------------------------------------------------------
        IF (xh(j).eq.0) then
! -----------------------------------------------------------
!	CORE: en_cond,  No anisotropy in the core
! -----------------------------------------------------------
          call en_cond(bm(i,j),tem,rho(j),ye(j),taun,taup,tcond,etab(i,j))
          f3=0.d0
          kt3=tcond*1.d13*1.d-40   ! cgs to [10^40 erg/(s*10^8 K*km)]
          kp3=kt3
          kh3=0.d0
        ELSE
! -----------------------------------------------------------
! 	impurity parameter (input for Potekhin's routine)
! -----------------------------------------------------------
          call get_Zimp (rho(j),Zimp)
! -----------------------------------------------------------
!	      CRUST: potekhin2019 (condconv, condegin)         
!    electron + phonons conductivity (not SF phonons)
! -----------------------------------------------------------
        call potekhin2019(tem,rho(j),bm(i,j),aa(j),zz(j),xh(j),Zimp,ye(j), &
     &  yn(j),taucru,tcond,tcondt,tcondh,rsigma,qj,qjt,qjh,CRKAPi,ksfph)

    
          if (bm(i,j).gt.0.d0) then
            f3=(tcond/tcondt-1.d0)/bm(i,j)**2
            kh3=(tcondh/tcondt)/bm(i,j)
          else
            f3=0.d0
            kh3=0.d0
          endif
          kt3=tcondt*1.d8*1.d5*1.d-40    ! cgs to [10^40 erg/(s*10^8 K*km)]
          etab(i,j)=CLIGHT**2/(4.d0*PI*rsigma)*T_YEAR*1.D-4   ! in Km^2/Myr
        ENDIF
!-----------------------------------------------------------------------
! Diffusion coefficients and anisotropy factors.
!-----------------------------------------------------------------------
        dfc2(k,l)=kt2
        dfc3(k,l)=kt3

        i=2*k-2
        j=2*l
        anis2(1,1,k,l)=1d0+br(i,j)**2*f2
        anis2(1,2,k,l)=bth(i,j)*br(i,j)*f2+bphi(i,j)*kh2
        anis2(2,1,k,l)=bth(i,j)*br(i,j)*f2-bphi(i,j)*kh2
        anis2(2,2,k,l)=1d0+bth(i,j)**2*f2

        i=2*k-1
        j=2*l-1
        anis3(1,1,k,l)=1d0+br(i,j)**2*f3
        anis3(1,2,k,l)=bth(i,j)*br(i,j)*f3+bphi(i,j)*kh3
        anis3(2,1,k,l)=bth(i,j)*br(i,j)*f3-bphi(i,j)*kh3
        anis3(2,2,k,l)=1d0+bth(i,j)**2*f3
       enddo
      enddo
!-----------------------------------------------------------------------
! Take averages of etab in the remaining locations (cell centers and corners)
! and fill the ghost cells just copying values of adjacent cells
! TODO: improve this, it's cumbersome
!-----------------------------------------------------------------------

      ! BC on ghost cells to be used in the loop
      etab(0,:) = etab(2,:)
      etab(nang+1,:) = etab(nang-1,:)
      do l=1,lmax
        do k=1,kmax
          j = 2*l
          i = 2*k-1
          etab(i,j)=0.5d0*(etab(i-1,j)+etab(i+1,j))
          if (k .ne. 1) then
            j = 2*l-1
            i = 2*k-2
            etab(i,j)=0.5d0*(etab(i-1,j)+etab(i+1,j))
          endif
        enddo
      enddo
      ! Repeat the BC on ghost cells since half of them have been filled
      etab(0,:) = etab(2,:)
      etab(nang+1,:) = etab(nang-1,:)
      etab(:,0) = etab(:,1)
      etab(:,np+1) = etab(:,np)
      etab(:,np+2) = etab(:,np)

    end subroutine compute_transport

end module microphysics
      
