!-----------------------------------------------------------------------
! Sample data for the stream functions P and T.
!
! Generate a textfile for chisq contours.
!-----------------------------------------------------------------------
	program SAMPLE

	implicit none

! Grid of the parameter space.
	integer i,j,ni,nj
	parameter(ni=100,nj=100)

! Fixed value of Pc.
	real*8 Pc
	parameter(Pc=0.35154804247367416d0)

! Poloidal and toroidal functions.
	integer n,ndim
	parameter(ndim=57)
	real*8, dimension(ndim) :: P,T,sig

! ma is the number of parameters.
! npc should be equal to or larger than ma (see LFIT or MRQMIN).
	integer ma,npc
	parameter(ma=3,npc=ma)
	integer ia(ma)
	real*8 a(ma)
	real*8 alamda,chisq
	real*8 alpha(npc,npc),covar(npc,npc)
	real*8 f,dfda(ma)
	external FTOR_NONL

! Number of iterations for the non-linear fit.
	integer niter
	parameter(niter=40)

! Read data.
	open(1,file='data.dat')
	do n=1,ndim
	   read(1,*)P(n),T(n)
	enddo
	close(1)
	sig=1d0

! Chisq map.
	open(2,file='chisq_map.txt')
	a(2)=Pc
	do i=1,ni
	   a(1)=2d0*dble(i-1)/dble(ni-1)
	   do j=1,nj
	      a(3)=0d0+2d0*dble(j-1)/dble(nj-1)
	      chisq=0d0
	      do n=1,ndim
		 call FTOR_NONL(P(n),a,f,dfda,ma)
		 chisq=chisq+(T(n)-f)**2
	      enddo
	      write(2,*)a(1),a(3),chisq
	   enddo
	   write(2,*)
	enddo
	close(2)

! The coefficients that should be kept fixed are indicated by a zero.
	ia(1)=1
	ia(2)=0
	ia(3)=1

! Initial guess.
	a(1)=1d0
	a(2)=Pc
	a(3)=1d0

! Non-linear fit.
	sig=1d0
	covar=0d0
	alamda=-1d0
	do i=1,niter
	    call MRQMIN(P,T,sig,ndim,a,ia,ma,covar,alpha,npc,chisq,
     &	    FTOR_NONL,alamda)
	enddo

	write(*,*)"Number of iterations (fixed):",niter
	write(*,*)"Chisq:",chisq
	write(*,*)"a(1):",a(1)
	write(*,*)"a(2):",a(2)
	write(*,*)"a(3):",a(3)

	stop
	end
!-----------------------------------------------------------------------
! Toroidal function that depends non-linearly on the parameters.
!-----------------------------------------------------------------------
	subroutine FTOR_NONL(x,a,f,dfda,ma)

	implicit none
	integer ma
	real*8 x,a(ma),f,dfda(ma)

! Function.
	f=a(1)*(x-a(2))**a(3)

! Derivatives wrt parameters.
	dfda(1)=(x-a(2))**a(3)
	dfda(2)=-a(1)*a(3)*(x-a(2))**(a(3)-1d0)
	dfda(3)=f*log(x-a(2))

	return
	end
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
	include "../nr_routines/covsrt.f"
	include "../nr_routines/gaussj.f"
	include "../nr_routines/lfit.f"
	include "../nr_routines/marquardt.f"
