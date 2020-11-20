!-----------------------------------------------------------------------
! Subroutine serving as an interface for ambipolar diffusion.
!
! Using B and J, calculate F, call AMBIPOLAR and return V.
! Called by BEVOL.
!-----------------------------------------------------------------------
      subroutine AMBIPOLAR_INTERFACE(lcore)
      use structure, only : fh
      use magnetic_evolution, only : br, bth, bphi, jr, jth, jphi
      implicit none
      include '../decl/dim2.h'


! Input variables.
! Needs dimensions kmax, lmax, lcore, np, nang.
! Needs br, bth, bphi, jr, jth, jphi.
! Needs fh.

! Output variables.
! Returns vr_amb, vth_amb, vphi_amb.

! Internally used variables.
	integer i,j,k,l,lcore
	real*8 fplus
	real*8 jrave,jthave

! Lorentz force.
	real*8, dimension (0:nang+1,0:np+2) :: fr,fth,fphi

! Variables needed by the ambipolar diffusion subroutine.
	real*8, dimension (3,lmax,kmax-1) :: fb_nc,v_amb
! Output variables
      real*8, dimension (0:nang+1,0:np+2) :: vr_amb,vth_amb,vphi_amb

!-----------------------------------------------------------------------
! Lorentz force density.
!-----------------------------------------------------------------------
	fr=0d0
	fth=0d0
	fphi=0d0

! Calculation of (odd,even) terms.
	do i=1,nang,2
	   do j=2,np,2
	      jrave=0.5d0*(jr(i,j+1)+jr(i,j-1))
	      jthave=0.5d0*(jth(i+1,j)+jth(i-1,j))

	      fr(i,j)=jthave*bphi(i,j)-jphi(i,j)*bth(i,j)
	      fth(i,j)=jphi(i,j)*br(i,j)-jrave*bphi(i,j)
	      fphi(i,j)=jrave*bth(i,j)-jthave*br(i,j)
	   enddo
	enddo

! Calculation of (even,odd) terms.
	do i=2,nang-1,2
	   do j=1,np-1,2
	      jrave=0.5d0*(jr(i+1,j)+jr(i-1,j))
	      jthave=0.5d0*(jth(i,j+1)+jth(i,j-1))

	      fr(i,j)=jthave*bphi(i,j)-jphi(i,j)*bth(i,j)
	      fth(i,j)=jphi(i,j)*br(i,j)-jrave*bphi(i,j)
	      fphi(i,j)=jrave*bth(i,j)-jthave*br(i,j)
	   enddo
	enddo

! Lorentz force per charge density.
	fb_nc=0d0
	do k=1,kmax-1
	   i=2*k
! Calculation at (even,odd) grid points.
	   do l=1,lcore
	      j=2*l-1
	      fB_nc(1,l,k)=fh(j)*fr(i,j)
	      fB_nc(2,l,k)=fh(j)*fth(i,j)
	      fB_nc(3,l,k)=fh(j)*fphi(i,j)
	   enddo
	enddo

!-----------------------------------------------------------------------
! Ambipolar diffusion velocity.
!-----------------------------------------------------------------------
	v_amb=0d0
	call AMBIPOLAR(lcore,fB_nc,v_amb)

! Amplification factor.
	fplus=1d2
	v_amb=v_amb*fplus

!-----------------------------------------------------------------------
! Output velocity.
!-----------------------------------------------------------------------
! Adapt the ambipolar diffusion velocity to the staggered grid.
	do k=1,kmax-1
	   i=2*k
! The (even,odd) terms are returned by AMBIPOLAR.
! Beyond the crust-core boundary the velocity is left as zero.
	   do l=1,lcore
	      j=2*l-1
	      vr_amb(i,j)=v_amb(1,l,k)
	      vth_amb(i,j)=v_amb(2,l,k)
	      vphi_amb(i,j)=v_amb(3,l,k)
	   enddo
	enddo

! Averages for the remaining terms.
	do i=2,nang-1,2
	   do j=1,2*lcore-1,2
! Averages for the (even,even) terms.
	      vr_amb(i,j+1)=(vr_amb(i,j+2)+vr_amb(i,j))/2d0
	      vth_amb(i,j+1)=(vth_amb(i,j+2)+vth_amb(i,j))/2d0
	      vphi_amb(i,j+1)=(vphi_amb(i,j+2)+vphi_amb(i,j))/2d0
! Averages for the (odd,odd) terms.
! (1,odd) terms are not calculated (and are left as zero).
	      vr_amb(i+1,j)=(vr_amb(i+2,j)+vr_amb(i,j))/2d0
	      vth_amb(i+1,j)=(vth_amb(i+2,j)+vth_amb(i,j))/2d0
	      vphi_amb(i+1,j)=(vphi_amb(i+2,j)+vphi_amb(i,j))/2d0
! Averages for the (odd,even) terms (using a four-way average).
! (1,even) terms are not calculated (and are left as zero).
	      vr_amb(i+1,j+1)=(vr_amb(i+2,j+2)+vr_amb(i+2,j)
     &	      +vr_amb(i,j+2)+vr_amb(i,j))/4d0
	      vth_amb(i+1,j+1)=(vth_amb(i+2,j+2)+vth_amb(i+2,j)
     &	      +vth_amb(i,j+2)+vth_amb(i,j))/4d0
	      vphi_amb(i+1,j+1)=(vphi_amb(i+2,j+2)+vphi_amb(i+2,j)
     &	      +vphi_amb(i,j+2)+vphi_amb(i,j))/4d0
	   enddo
	enddo

! Manually set values at the boundaries.
! Along the north axis.
	vr_amb(1,:)=vr_amb(2,:)
	vth_amb(1,:)=0d0
	vphi_amb(2,:)=0d0
! Along the south axis.
	vr_amb(nang,:)=vr_amb(nang-1,:)
	vth_amb(nang,:)=0d0
	vphi_amb(nang,:)=0d0
! At the crust-core boundary.
	vr_amb(:,2*lcore)=0d0
	vth_amb(:,2*lcore)=vth_amb(:,2*lcore-1)
	vphi_amb(:,2*lcore)=vphi_amb(:,2*lcore-1)

	return
	end
