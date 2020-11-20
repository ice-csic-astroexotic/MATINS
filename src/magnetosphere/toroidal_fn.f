!-----------------------------------------------------------------------
! Toroidal functions for the force-free magnetosphere.
! Needed by FORCE_FREE, INTERFACE and TWIST.
!
! Contents:
! FTOR
! FTOR_LIN
! DFTOR_LIN
! FTOR_NONL
! FTOR_NONL_SIGMA1
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Evaluate the toroidal function and its derivative.
!-----------------------------------------------------------------------
	subroutine FTOR(x,a,f,dfdx,ma,fit_type)

	implicit none
	integer fit_type,ma
	real*8 x,f,dfdx
	real*8 a(ma),fn(ma),dfn(ma)

! Linear fit for the toroidal function.
	if(fit_type.eq.1)then
	   call FTOR_LIN(x,fn,ma)
	   call DFTOR_LIN(x,dfn,ma)
	   f=DOT_PRODUCT(a,fn)
	   dfdx=DOT_PRODUCT(a,dfn)

! Non-linear fit for the toroidal function.
	else if(fit_type.eq.2)then
	   call FTOR_NONL(x,a,f,dfn,ma)
! Note that dfdx=-dfn(2) when f is a function of x-a(2).
	   dfdx=-dfn(2)
	endif

	return
	end
!-----------------------------------------------------------------------
! Toroidal function that depends linearly on the parameters.
! The corresponding derivatives must the correctly defined in DFTOR_LIN.
!-----------------------------------------------------------------------
	subroutine FTOR_LIN(x,f,ma)

	implicit none
	integer i,ma
	real*8 x,f(ma)

! Variables passed on from INTERFACE.
	integer fn_type
	common /function_type/ fn_type
	real*8 Pc
	common /p_critical/ Pc

	if(fn_type.eq.1)then
! Polynomial functions, not guaranteed to give T=0 at P=Pc.
	   do i=1,ma
	      f(i)=x**i
	   enddo
	else if(fn_type.eq.2)then
! Power-law functions, guaranteed to give T=0 at P=Pc.
	   do i=1,ma
	      f(i)=(x-Pc)**i
	   enddo
	else
	   write(*,'(a)')"FTOR_LIN: Unexpected value of fn_type."
	   stop
	endif

	return
	end
!-----------------------------------------------------------------------
! Derivatives of the fitting functions defined in FTOR_LIN.
!-----------------------------------------------------------------------
	subroutine DFTOR_LIN(x,df,ma)

	implicit none
	integer i,ma
	real*8 x,df(ma)

! Variables passed on from INTERFACE.
	integer fn_type
	common /function_type/ fn_type
	real*8 Pc
	common /p_critical/ Pc

	if(fn_type.eq.1)then
! Polynomial functions, not guaranteed to give T=0 at P=Pc.
	   do i=1,ma
	      df(i)=dble(i)*x**(i-1)
	   enddo
	else if(fn_type.eq.2)then
! Power-law functions, guaranteed to give T=0 at P=Pc.
	   do i=1,ma
	      df(i)=dble(i)*(x-Pc)**(i-1)
	   enddo
	endif

	return
	end
!-----------------------------------------------------------------------
! Toroidal function that depends non-linearly on the parameters.
!-----------------------------------------------------------------------
	subroutine FTOR_NONL(x,a,f,dfda,ma)

	implicit none
	integer ma
	real*8 x,a(ma),f,dfda(ma)

! Variables passed on from INTERFACE.
	integer fn_type
	common /function_type/ fn_type

! a(2)=Pc.
	if(x.gt.a(2))then
	   if(fn_type.eq.1)then
! Function.
	      f=a(1)*(x-a(2))**a(3)
! Derivatives wrt parameters.
	      dfda(1)=(x-a(2))**a(3)
	      dfda(2)=-a(1)*a(3)*(x-a(2))**(a(3)-1d0)
	      dfda(3)=f*log(x-a(2))
	   else if(fn_type.eq.2)then
! Function.
	      f=a(1)*(x-a(2))+a(3)*(x-a(2))**2
! Derivatives wrt parameters.
	      dfda(1)=x-a(2)
	      dfda(2)=-a(1)-2d0*a(3)*(x-a(2))
	      dfda(3)=(x-a(2))**2
	   else
	      write(*,'(a)')"FTOR_NONL: Unexpected value of fn_type."
	      stop
	   endif
	else
	   f=0d0
	   dfda=0d0
	endif

	return
	end
!-----------------------------------------------------------------------
! Special case of the non-linear function for fn_type=1.
!-----------------------------------------------------------------------
	subroutine FTOR_NONL_SIGMA1(x,f,ma)

	implicit none
	integer ma
	real*8 x,f(ma)

! Variables passed on from INTERFACE.
	real*8 Pc
	common /p_critical/ Pc

! Function.
	f(1)=x-Pc
	f(2)=1d0
	f(3)=1d0

	return
	end
