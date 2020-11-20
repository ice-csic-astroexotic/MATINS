!-------------------------------------------------------------------------------
! Magneto Thermal 2D
!-------------------------------------------------------------------------------
! 
!
!-------------------------------------------------------------------------------


!-----------------------------------------------------------------------
! Contents:
! RECONSTRUCTION (used in 2012 code)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! New version of the subroutine.
!-----------------------------------------------------------------------
  subroutine RECONSTRUCTION(xc,xp,xm,xright,xleft,
     &	yc,yp,ym,yright,yleft,recon)

	implicit none
	real*8 xc,xp,xm,yc,yp,ym
	real*8 xright,xleft,xhalf
	real*8 alpham,alphap,alphac,alpha
	real*8 yright,yleft
	integer recon

	alphap=(yp-yc)/(xp-xc)
	alpham=(yc-ym)/(xc-xm)
	alphac=(yp-ym)/(xp-xm)
	if(alpham*alphap.le.0d0)then
	   alpha=0d0
	else
! Method: minmod (recon=1) or mc (recon=2).
	   if(recon.eq.1)then
	      if(alpham.gt.0d0)then
		 alpha=min(alpham,alphap)
	      else
		 alpha=max(alpham,alphap)
	      endif
	   else
	      if(alpham.gt.0d0)then
		 alpha=min(2d0*alpham,alphac,2d0*alphap)
	      else
		 alpha=max(2d0*alpham,alphac,2d0*alphap)
	      endif
	   endif
	endif

! Old version.
cc	xhalf=0.5d0*(xc+xp)
cc	yright=yc+alpha*(xhalf-xc)
cc	xhalf=0.5d0*(xc+xm)
cc	yleft=yc+alpha*(xhalf-xc)

! New version.
	yright=yc+alpha*(xright-xc)
	yleft=yc+alpha*(xleft-xc)

	return
	end
