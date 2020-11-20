!-----------------------------------------------------------------------
! Fluxtube velocity in a type II superconducting core.
! (Following the notes of Glampedakis 2012.)
!
! Inverse curvature radius:
! k = (b dot grad) b, b unit vector, so k = - b x curl(b).
! velocity in km/Myr, fluxen=1d7 erg/cm, quantum=2d-3 cm2/s.
! Rp defined in input.f, together with rhop.
!-----------------------------------------------------------------------
	subroutine BFLUXTUBE

	implicit none
	include '../decl/dim2.h'
	include '../decl/varsB.h'
	include '../decl/adv.h'

! Output variables: vrexp,vthexp,vphiexp.

! Internally used variables.
	integer i,j
	real*8 a,bvert,fluxen,kmmyr,quantum
	real*8 curlr,curlth,curlphi,fr,fth,fphi
	parameter(fluxen=1d7,kmmyr=1d6*yrs/1d5,quantum=2d-3)

	do j=2,2*lcore
	   a=fluxen*rp(j)*kmmyr/(quantum*rhop(j))/(1+rp(j)**2)
	   do i=2,nang-1
	      curlr=(bphin(i+1,j)*sth(i+1)
     &	      -bphin(i-1,j)*sth(i-1))/(lth(j)*sth(i))
	      curlth=-(bphin(i,j+1)*rb(j+1)
     &	      -bphin(i,j-1)*rb(j-1))/(lr(j)*rb(j))
	      curlphi=(rb(j+1)*bthn(i,j+1)
     &	      -rb(j-1)*bthn(i,j-1))/(lr(j)*rb(j))
     &	      -(brn(i+1,j)-brn(i-1,j))/lth(j)

! Lorentz force.
	      fr=curlth*bphin(i,j)-curlphi*bthn(i,j)
	      fth=curlphi*brn(i,j)-curlr*bphin(i,j)
	      fphi=curlr*bthn(i,j)-curlth*brn(i,j)

! Buoyancy (just the leading term included).
	      fr=fr+(1d0-belam(j)**(-2))/(2d0*rb(j)*cs2(j/2))

! Velocity field.
	      vrexp(i,j)=a*(-rp(j)*curlr+fr)
	      vthexp(i,j)=a*(-rp(j)*curlth+fth)
	      vphiexp(i,j)=a*(-rp(j)*curlphi+fphi)
	   enddo
	enddo

	return
	end
