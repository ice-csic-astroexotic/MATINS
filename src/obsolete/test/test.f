!-----------------------------------------------------------------------
! Test for linear and non-linear least squares fit with LFIT and MRQMIN.
!-----------------------------------------------------------------------
	program TEST_MRQMIN

	implicit none

! Number of parameters and data points.
	integer ma,nca,ndata
	parameter(ma=3,nca=3,ndata=10)

! Variables passed on to MRQMIN.
	integer ia(ma)
	real*8 alamda,chisq,a(ma),alin(ma)
	real*8 alpha(nca,nca),covar(nca,nca)
	real*8 sig(ndata),x(ndata),y(ndata)
	external FGAUSS,FLINEAR,FPOL,FTOR

! Other internally used variables.
	integer i
	real*8 pi
	real*8 Pc,s,sigma
	real*8 f,flin,dfda(ma)

! Number of iterations.
	integer iter
	parameter(iter=20)

! Definition of pi.
	pi=2d0*asin(1d0)

! Data.
	Pc=0.5d0
	s=2d0
	sigma=2d0
	do i=1,ndata
	   x(i)=Pc+(1d0-Pc)*dble(i)/dble(ndata)
	   y(i)=s*(x(i)-Pc)**sigma
	enddo
	sig=1d0

! Parameters to be fitted for (0 in order to keep fixed).
	ia=1

! Initial guess for the parameters.
	a=0.5d0
	alin=a

! Linear least squares fit.
	call LFIT(x,y,sig,ndata,alin,ia,ma,covar,nca,chisq,FLINEAR)

! Non-linear least squares fit.
	alamda=-1d0
	do i=1,iter
	   call MRQMIN(x,y,sig,ndata,a,ia,ma,covar,alpha,nca,chisq,
     &	   FTOR,alamda)
	enddo

! Output.
	do i=1,ndata
	   call FLINEAR(x(i),dfda,ma)
	   flin=alin(1)*dfda(1)+alin(2)*dfda(2)+alin(3)*dfda(3)
	   call FPOL(x(i),a,f,dfda,ma)
	   write(1,*)x(i),y(i),f,flin
	enddo

	stop
	end
!-----------------------------------------------------------------------
! Gaussian functions given as an example in Numerial Recipes.
!-----------------------------------------------------------------------
	subroutine FGAUSS(x,a,y,dyda,na)

	implicit none
	integer na
	real*8 x,y,a(na),dyda(na)

	integer i
	real*8 arg,ex,fac

	y=0d0
	do i=1,na-1,3
	   arg=(x-a(i+1))/a(i+2)
	   ex=exp(-arg**2)
	   fac=a(i)*ex*2d0*arg
	   y=y+a(i)*ex
	   dyda(i)=ex
	   dyda(i+1)=fac/a(i+2)
	   dyda(i+2)=fac*arg/a(i+2)
	enddo

	return
	end
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
	subroutine FLINEAR(x,dyda,na)

	implicit none
	integer na
	real*8 x,dyda(na)

	dyda(1)=x
	dyda(2)=x**2
	dyda(3)=x**3

	return
	end
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
	subroutine FPOL(x,a,y,dyda,na)

	implicit none
	integer na
	real*8 x,y,a(na),dyda(na)

	y=a(1)*x+a(2)*x**2+a(3)*x**3
	dyda(1)=x
	dyda(2)=x**2
	dyda(3)=x**3

	return
	end
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
	subroutine FTOR(x,a,y,dyda,na)

	implicit none
	integer na
	real*8 x,y,a(na),dyda(na)

	y=a(1)*(x-a(2))**a(3)
	dyda(1)=(x-a(2))**a(3)
	dyda(2)=-a(1)*a(3)*(x-a(2))**(a(3)-1d0)
	dyda(3)=y*log(a(1)*(x-a(2)))

	return
	end
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
	include "../nr_routines/covsrt.f"
	include "../nr_routines/gaussj.f"
	include "../nr_routines/lfit.f"
	include "../nr_routines/marquardt.f"
