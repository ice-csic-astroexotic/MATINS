!!-------------------------------------------------------------------------------
!! Module: Microphysics
!!
!!> @brief This module is responsible for all the microphysics input.
!!
!!-------------------------------------------------------------------------------
module conductivities

  use constants, only: PI, CLIGHT, T_YEAR, UNIT_T, UNIT_R, UNIT_EN, UNIT_TIME
  use grid, only: nangt, nrt, nang, nr, ncore, r
  use grid, only: tem0, temp, bm, etab, omegatau_arr, kappa_perp_arr
  use grid, only: rho, ne, nn, npr, kFe, kFn, kFp, effme, effmn, effmp
  use grid, only: ye, yn, xh, aa, zz, tccru, Zimp

  contains

!!-----------------------------------------------------------------------
!!   Interface with the microphysics routines
!!   Calculates diffusion coefficients, anisotropy factors and resistivity.
!!
!! Contents:
!! compute_conductivities
!!-----------------------------------------------------------------------

    subroutine compute_conductivities

      ! Modules

      implicit none

      ! Internally used variables.
      integer i,j,k,p,it,jt,kt,is
      real*8 ff,kh
      real*8 tem, taucru
      real*8 tcond, tcondperp, tcondhall, rsigma

      kappa_perp_arr = 0d0
      omegatau_arr = 0d0
      etab = 0d0

      do it=1,nrt   ! Start loop from CC interface
       i = it*2 - 1  ! odd values
       is = i + ncore
       do p=1,6  
        do jt=0,nangt+1
          j = jt*2
          do kt=0,nangt+1
            k = kt*2

            tem = tem0(it,jt,kt,p)
            ! T/Tc for SF, with a floor to avoid divisions by 0 or small numbers
            taucru = tem/dmax1(tccru(is),0.1d0*tem)
            
            ! crust: electron + phonons conductivity
            call cond_crust(tem,rho(is),bm(i,j,k,p),aa(is),zz(is),xh(is),Zimp(i),ye(is), &
            &                 yn(is),kFn(is),effmn(is),taucru,tcond,tcondperp,tcondhall, &
            &                 rsigma)

            ! Magnetic diffusivity, from rsigma, which is the electrical conductivity
            etab(i,j,k,p)=CLIGHT**2/(4.d0*PI*rsigma)*T_YEAR*1.D-4   ! in Km^2/Myr

            ! Here we use the notation of e.g. Aguilera et al. 2008
            ! However, we approximate with the classical formulation,
            ! hence the perpendicular conductivity and omegatau are enough,
            ! no need for the Hall conductivity
            if (bm(i,j,k,p) > 0.d0) then
              ff = (tcond/tcondperp-1.d0)/bm(i,j,k,p)**2 ! ff = (omegatau/B)^2
            !  kh = (tcondhall/tcondperp)/bm(i,j,k,p)
            else
              ff = 0.d0
            !  kh = 0.d0
            endif
            kappa_perp_arr(it,jt,kt,p) = tcondperp*UNIT_TIME*UNIT_T*UNIT_R/UNIT_EN  ! cgs to [10^40 erg/(Myr*10^8 K*km)]
            omegatau_arr(it,jt,kt,p) = dsqrt(ff)*bm(i,j,k,p)

          enddo
        enddo
       enddo
      enddo
    
      ! Take averages of etab at all the other places
      ! Interface centers
      etab(2:nr-2:2,0:nang+1:2,0:nang+1:2,:) = 0.5*(etab(1:nr-3:2,0:nang+1:2,0:nang+1:2,:) + etab(3:nr-1:2,0:nang+1:2,0:nang+1:2,:))  ! Center radial interfaces
      etab(:,1:nang:2,0:nang+1:2,:) = 0.5*(etab(:,0:nang-1:2,0:nang+1:2,:) + etab(:,2:nang+1:2,0:nang+1:2,:))  ! Center xi interfaces
      etab(:,0:nang+1:2,1:nang:2,:) = 0.5*(etab(:,0:nang+1:2,0:nang-1:2,:) + etab(:,0:nang+1:2,2:nang+1:2,:))  ! Center eta interfaces

      ! Edges middle points (using the average among two interface centers, it's the fastest way)
      etab(2:nr-2:2,3:nang-2:2,2:nang-1:2,:) = 0.5*(etab(1:nr-3:2,3:nang-2:2,2:nang-1:2,:) + etab(3:nr-1:2,3:nang-2:2,2:nang-1:2,:))  ! eta-direction edge
      etab(2:nr-2:2,2:nang-1:2,3:nang-2:2,:) = 0.5*(etab(1:nr-3:2,2:nang-1:2,3:nang-2:2,:) + etab(3:nr-1:2,2:nang-1:2,3:nang-2:2,:))  ! xi-direction edge
      etab(1:nr-1:2,3:nang-2:2,3:nang-2:2,:) = 0.5*(etab(1:nr-1:2,2:nang-3:2,3:nang-2:2,:) + etab(1:nr-1:2,4:nang-1:2,3:nang-2:2,:))  ! r-direction edge
      
      ! Vertexes (using the average among two edges middle points)
      etab(2:nr-2:2,3:nang-2:2,3:nang-2:2,:) = 0.5*(etab(1:nr-3:2,3:nang-2:2,3:nang-2:2,:) + etab(3:nr-1:2,3:nang-2:2,3:nang-2:2,:))

      ! Special treatment between patches
      ! Equatorial edges
      ! patches I and II, II and III and III and IV
      do p=1,3 
        etab(:,nang,:,p) = 0.5*(etab(:,nang-1,:,p) + etab(:,2,:,p+1))
        etab(:,1,:,p+1) = etab(:,nang,:,p)
      enddo
      ! patches IV and I
      etab(:,nang,:,4) = 0.5*(etab(:,nang-1,:,4) + etab(:,2,:,1))
      etab(:,1,:,1) = etab(:,nang,:,4)
      
      ! edges - to be double checked later on- nothing to do with ghost cells so far 

      ! North edges
      ! patches I and V
      etab(:,:,nang,1) = 0.5*(etab(:,:,nang-1,1) + etab(:,:,2,5))
      etab(:,:,1,5) = etab(:,:,nang,1)
      ! patches II and V
      etab(:,:,nang,2) = 0.5*(etab(:,:,nang-1,2) + etab(:,nang-1,:,5))
      etab(:,nang,:,5) = etab(:,:,nang,2)      
      ! patches III and V
      etab(:,1:nang,nang,3) = 0.5*(etab(:,1:nang,nang-1,3) + etab(:,nang:1:-1,nang-1,5))
      etab(:,nang:1:-1,nang,5) = etab(:,1:nang,nang,3)     
      ! patches IV and V
      etab(:,1:nang,nang,4) = 0.5*(etab(:,1:nang,nang-1,4) + etab(:,2,nang:1:-1,5))
      etab(:,1,nang:1:-1,5) = etab(:,1:nang,nang,4)      

      ! South edges
      ! patches I and VI
      etab(:,:,1,1) = 0.5*(etab(:,:,2,1) + etab(:,:,nang-1,6))
      etab(:,:,nang,6) = etab(:,:,1,1)
      ! patches II and VI
      etab(:,1:nang,1,2) = 0.5*(etab(:,1:nang,2,2) + etab(:,nang-1,nang:1:-1,6))
      etab(:,nang,nang:1:-1,6) = etab(:,1:nang,1,2)
      ! patches III and VI
      etab(:,1:nang,1,3) = 0.5*(etab(:,1:nang,2,3) + etab(:,nang:1:-1,2,6))
      etab(:,nang:1:-1,1,6) = etab(:,1:nang,1,3)   
      ! patches IV and VI
      etab(:,:,1,4) = 0.5*(etab(:,:,2,4) + etab(:,2,:,6))
      etab(:,1,:,6) = etab(:,:,1,4)
      
      ! Vertices 
      ! Among patches 1,4,6
      etab(:,1,1,1) = 1d0/3d0*(etab(:,2,1,1) + etab(:,1,2,1) + etab(:,nang-1,1,4))
      etab(:,nang,1,4) = etab(:,1,1,1)
      etab(:,1,nang,6) = etab(:,1,1,1)
      ! Among patches 2,1,6
      etab(:,1,1,2) = 1d0/3d0*(etab(:,2,1,2) + etab(:,1,2,2) + etab(:,nang-1,1,1))
      etab(:,nang,1,1) = etab(:,1,1,2)
      etab(:,nang,nang,6) = etab(:,1,1,2)
      ! Among patches 3,2,6
      etab(:,1,1,3) = 1d0/3d0*(etab(:,2,1,3) + etab(:,1,2,3) + etab(:,nang-1,1,2))
      etab(:,nang,1,2) = etab(:,1,1,3)
      etab(:,nang,1,6) = etab(:,1,1,3)
      ! Among patches 4,3,6
      etab(:,1,1,4) = 1d0/3d0*(etab(:,2,1,4) + etab(:,1,2,4) + etab(:,nang-1,1,3))
      etab(:,nang,1,3) = etab(:,1,1,4)
      etab(:,1,1,6) = etab(:,1,1,4)
      ! Among patches 1,4,5
      etab(:,1,nang,1) = 1d0/3d0*(etab(:,2,nang,1) + etab(:,1,nang-1,1) + etab(:,nang-1,nang,4))
      etab(:,nang,nang,4) = etab(:,1,nang,1)
      etab(:,1,1,5) = etab(:,1,nang,1)
      ! Among patches 2,1,5
      etab(:,1,nang,2) = 1d0/3d0*(etab(:,2,nang,2) + etab(:,1,nang-1,2) + etab(:,nang-1,nang,1))
      etab(:,nang,nang,1) = etab(:,1,nang,2)
      etab(:,nang,1,5) = etab(:,1,nang,2)
      ! Among patches 3,2,5
      etab(:,1,nang,3) = 1d0/3d0*(etab(:,2,nang,3) + etab(:,1,nang-1,3) + etab(:,nang-1,nang,2))
      etab(:,nang,nang,2) = etab(:,1,nang,3)
      etab(:,nang,nang,5) = etab(:,1,nang,3)
      ! Among patches 4,3,5
      etab(:,1,nang,4) = 1d0/3d0*(etab(:,2,nang,4) + etab(:,1,nang-1,4) + etab(:,nang-1,nang,3))
      etab(:,nang,nang,3) = etab(:,1,nang,4)
      etab(:,1,nang,5) = etab(:,nang,1,4)

    ! Radial BC, just copying the last available value
      etab(nr,:,:,:) = etab(nr-1,:,:,:) !+ (r(nr)-r(nr-1))*(etab(nr-1,:,:,:)-etab(nr-2,:,:,:)) &
    ! & /(r(nr-1)-r(nr-2))

      etab(0,:,:,:) = etab(1,:,:,:)  ! This is not really needed
      etab(nr+1,:,:,:) = etab(nr,:,:,:)  ! This is not really needed

    end subroutine compute_conductivities

end module conductivities
