!-----------------------------------------------------------------------
! Numerical Recipes subroutine for associated Legendre polynomials.
!
! Contents:
! PLGNDR
! DPLGNDR
!-----------------------------------------------------------------------
	double precision function PLGNDR(l,m,x)

	implicit none
	integer l,m
	double precision x
	integer i,ll
	double precision fact,pll,pmm,pmmp1,somx2

	if(m.lt.0.or.m.gt.l.or.dabs(x).gt.1d0)then
	   stop "PLGNDR: Bad arguments!"
	endif

	pmm=1d0
	if(m.gt.0)then
	   somx2=dsqrt((1d0-x)*(1d0+x))
	   fact=1d0
	   do i=1,m
	      pmm=-pmm*fact*somx2
	      fact=fact+2d0
	   enddo
	endif

	if(l.eq.m)then
	   PLGNDR=pmm
	else
	   pmmp1=x*dble(2*m+1)*pmm
	   if(l.eq.m+1)then
	      PLGNDR=pmmp1
	   else
	      do ll=m+2,l
		 pll=(x*dble(2*ll-1)*pmmp1-dble(ll+m-1)*pmm)/dble(ll-m)
		 pmm=pmmp1
		 pmmp1=pll
	      enddo
	      PLGNDR=pll
	   endif
	endif

	return
	end
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
	double precision function DPLGNDR(l,m,x)

	implicit none
	integer l,m
	double precision x,PLGNDR

	if(dabs(x).eq.1d0)then
	   stop "DPLGNDR: Derivatives at x=+/-1 not implemented!"
	endif

	DPLGNDR=(dble(l)*x*PLGNDR(l,m,x)-dble(l+m)*PLGNDR(l-1,m,x))/
     &	(x**2d0-1d0)

	return
	end
