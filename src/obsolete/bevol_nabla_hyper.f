
!-----------------------------------------------------------------------
! Hyper-resistivity term used by BEVOL.
!-----------------------------------------------------------------------
	subroutine NABLA_HYPER(f,nabla4)

	implicit none
	include '../decl/dim2.h'

! Input and output variables.
	real*8 f(0:kd+1,0:ltot+2),nabla4(kd,ltot+1)

! Internally used variables.
	integer k,l
	real*8 f03,f04,f10,f12,f20,f22,f30,f40
	real*8 csc,cotg,dr(ltot),dz

! Boundary values for the input function f.
	f(:,ltot+1)=0d0
	f(:,ltot+2)=0d0
	f(:,0)=0d0
	f(:,1)=0d0
	f(1,:)=-f(2,:)
	f(0,:)=-f(3,:)
	f(kd,:)=-f(kmax,:)
	f(kd+1,:)=-f(kmax-1,:)

! Differential elements of the grid.
	dz=zc(3)-zc(2)
	do l=2,ltot
	   dr(l)=0.5d0*(r_ext(l+1)-r_ext(l-1))
	enddo
	dr(1)=dr(2)

! Calculation of nabla4.
! Warning: This calculation is valid for an equispaced grid!
! Warning: The radii do not correspond to the points where f is given.
	do k=2,kmax
	   csc=(1d0/dsin(zc(k)))**2
	   cotg=dcos(zc(k))/dsin(zc(k))
	   do l=lc+2,ltot

	      f10=0.5d0*(f(k+1,l)-f(k-1,l))/dz
	      f20=(f(k+1,l) - 2d0*f(k,l) + f(k-1,l) )/dz**2
	      f30=(0.5d0*(f(k+2,l)-f(k-2,l))-f(k+1,l)+f(k-1,l))/dz**3
	      f40=(f(k+2,l)+f(k-2,l)-4d0*(f(k+1,l)+f(k-1,l))
     &		 +6d0*f(k,l))/dz**4
	      f03=(0.5d0*(f(k,l+2)-f(k,l-2))-f(k,l+1)+f(k,l-1))/dr(l)**3
	      f04=(f(k,l+2)+f(k,l-2)-4d0*(f(k,l+1)+f(k,l-1))
     &		 +6d0*f(k,l))/dr(l)**4

	      f12 = 0.5d0/(dz*dr(l)**2)*
     &		 ( ( f(k+1,l+1) - 2d0*f(k+1,l) + f(k+1,l-1))
     &		 - ( f(k-1,l+1) - 2d0*f(k-1,l) + f(k-1,l-1)) )

	      f22 = 1d0/(dz**2*dr(l)**2)*
     &		 ( f(k+1,l+1) - 2d0*f(k+1,l) + f(k+1,l-1)
     &		 + f(k-1,l+1) - 2d0*f(k-1,l) + f(k-1,l-1)
     &		 - 2d0*f(k,l+1) + 4d0*f(k,l) - 2d0*f(k,l-1) )

	      nabla4(k,l) =
     &		 1d0/r_ext(l)**4*(cotg*(2d0+csc)*f10
     &		 - cotg**2*f20 + 2d0*cotg*f30 + f40)
     &		 + 1d0/r_ext(l)**2*(2d0*cotg*f12 + 2d0*f22)
     &		 + 4d0/r_ext(l)*f03 + f04

	      nabla4(k,l)=nabla4(k,l)
     &		 /(1d0/(r_ext(l)*dz)**2+1d0/dr(l)**2)
	   enddo
	enddo

! Boundary values for nabla4.
	nabla4(:,lc+1)=0d0
	nabla4(:,lc)=0d0

	return
	end
!-----------------------------------------------------------------------
! New routine for hyper-resistivity.
!-----------------------------------------------------------------------
	subroutine NABLA_HYPER_NEW(kd,ld,zc,rc,f,nabla4)

	implicit none
	integer k,l,kd,ld
	real*8, dimension (kd) :: zc
	real*8, dimension (ld) :: rc
	real*8, dimension (kd,ld) :: f,nabla2,nabla4
	real*8 dr,dz

	call LAPLACIAN(kd,ld,zc,rc,f,nabla2)
	call LAPLACIAN(kd,ld,zc,rc,nabla2,nabla4)

! Rescale by grid size.
	do k=2,kd
	   dz=zc(k)-zc(k-1)
	   do l=2,ld
	      dr=rc(l)-rc(l-1)
	      nabla4(k,l)=nabla4(k,l)/(1d0/dr**2+1d0/(rc(l)*dz)**2)
	   enddo
! Rescale the boundary for l=1.
	   dr=rc(2)-rc(1)
	   nabla4(k,1)=nabla4(k,1)/(1d0/dr**2+1d0/(rc(1)*dz)**2)
	enddo
! Rescale the boundary for k=1 and l=1.
	dz=zc(2)-zc(1)
	nabla4(1,1)=nabla4(1,1)/(1d0/dr**2+1d0/(rc(1)*dz)**2)
! Rescale the boundary for k=1.
	do l=2,ld
	   dr=rc(l)-rc(l-1)
	   nabla4(1,l)=nabla4(1,l)/(1d0/dr**2+1d0/(rc(l)*dz)**2)
	enddo

	return
	end
!-----------------------------------------------------------------------
! Laplacian of a function.
!-----------------------------------------------------------------------
	subroutine LAPLACIAN(kd,ld,zc,rc,f,nabla2)

	implicit none
	integer k,l,kd,ld
	real*8, dimension (kd,ld) :: f,nabla2
	real*8, dimension (ld) :: rc,a1,b1,c1,a2,b2,c2
	real*8, dimension (kd) :: zc,za1,zb1,zc1,za2,zb2,zc2
	real*8 alpha,dfr1,dfr2,dfc1,dfc2, cotgz

	nabla2=0d0

! Coefficients used for the derivatives in uneven grids.
	call COEFFD2(ld,rc,a1,b1,c1,a2,b2,c2)
	call COEFFD2(kd,zc,za1,zb1,zc1,za2,zb2,zc2)

	do k=2,kd-1
	   do l=2,ld-1
	      dfr1 = a1(l)*f(k,l-1)+b1(l)*f(k,l)+c1(l)*f(k,l+1)
	      dfr2 = a2(l)*f(k,l-1)+b2(l)*f(k,l)+c2(l)*f(k,l+1)

	      dfc1 = za1(k)*f(k-1,l)+zb1(k)*f(k,l)+zc1(k)*f(k+1,l)
	      dfc2 = za2(k)*f(k-1,l)+zb2(k)*f(k,l)+zc2(k)*f(k+1,l)

	      cotgz=dcos(zc(k))/dsin(zc(k))
	      nabla2(k,l) = dfr2 + 2d0*dfr1/rc(l)
     &	      + (dfc2+cotgz*dfc1)/rc(l)**2
! The extra term that appears in the vector Laplacian (vs. the scalar).
! Seems to have a negligible effect!
     &	      - f(k,l)/(rc(l)*sin(zc(k)))**2
	   enddo
	enddo

! Linear extrapolation for the radial boundaries.
	alpha=(rc(2)-rc(1))/(rc(3)-rc(2))
	nabla2(:,1)=nabla2(:,2)*(1d0+alpha)-alpha*nabla2(:,3)

	alpha=(rc(ld)-rc(ld-1))/(rc(ld-1)-rc(ld-2))
	nabla2(:,ld)=nabla2(:,ld-1)*(1d0+alpha)-alpha*nabla2(:,ld-2)

! Conditions at the axis.
	nabla2(1,:)=0.d0
	nabla2(kd,:)=0.d0

	return
	end
!-----------------------------------------------------------------------
! Coefficients used for the derivatives in uneven grids.
!-----------------------------------------------------------------------
	subroutine COEFFD2(ld,r,a1,b1,c1,a2,b2,c2)

	implicit none
	integer l,ld
	real*8, dimension (ld) :: r,a1,b1,c1,a2,b2,c2

	a1=0d0
	b1=0d0
	c1=0d0
	a2=0d0
	b2=0d0
	c2=0d0
	do l=2,ld-1
	   a1(l)=-(r(l+1)-r(l))/((r(l)-r(l-1))*(r(l+1)-r(l-1)))
	   b1(l)=(r(l+1)-2d0*r(l)+r(l-1))/((r(l)-r(l-1))*(r(l+1)-r(l)))
	   c1(l)=(r(l)-r(l-1))/((r(l+1)-r(l))*(r(l+1)-r(l-1)))
	   a2(l)=2d0/((r(l)-r(l-1))*(r(l+1)-r(l-1)))
	   b2(l)=-2d0/((r(l)-r(l-1))*(r(l+1)-r(l)))
	   c2(l)=2d0/((r(l+1)-r(l))*(r(l+1)-r(l-1)))
	enddo

	return
	end

!-----------------------------------------------------------------------
! Mean square averages and unit vectors of magnetic field components.
! Completing the grid of components where not reconstructed.
!-----------------------------------------------------------------------
	subroutine BMODULUS

	implicit none
	include '../decl/dim2.h'
	include '../decl/varsB.h'
	integer i,j,k,l

! Magnitude of the magnetic field.
        bm=dsqrt(br**2+bth**2+bphi**2)

! Magnitude in the reduced (original) grid.
	bmed=0d0
	do k=1,kmax
	   i=2*k-2
	   do l=lc,lmax
	      j=2*l-1
	      bmed(k,l)=dsqrt(br(2*k-2,2*l-1)**2+bth(2*k-2,2*l-1)**2+bphi(2*k-2,2*l-1)**2)
	   enddo
	enddo

	return
	end

