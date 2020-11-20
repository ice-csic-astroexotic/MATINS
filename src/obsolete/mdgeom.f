!-----------------------------------------------------------------------
! Calculates geometrical arrays for the diffusion scheme.
!
! Definitions of variables:
! r	Radial coordinate of Lagrangian-mesh corners.
! z	Angular coordinate of Lagrangian-mesh corners.
! kmax	Maximum angular index (k) of real cells and edges.
! lmax	Maximum radial index (l) of real cells and edges.
! kd	k-dimensioning parameter.
! ld	l-dimensioning parameter.
!-----------------------------------------------------------------------
	subroutine MDGEOM()
	
	use constants, only: PI
	implicit none
	include '../decl/dim2.h'
	integer k,l


! Angular grids.
! z(kd) a ghost cell, and is probably not used!
	z=0d0
	do k=1,kd
	   z(k)=pi*dble(k-1)/dble(kmax-1)
	enddo

c	zc=0d0
c	do k=2,kd
c	   zc(k)=0.5d0*(z(k)+z(k-1))
c	enddo
c	zc(1)=-zc(2)

! Centered radial grid.
c	rc=0d0
c	rc(1)=0.5d0*r(1)
c	do l=2,ld
c	   rc(l)=0.5d0*(r(l)+r(l-1))
c	enddo

	return
	end
