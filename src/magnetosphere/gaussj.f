!-----------------------------------------------------------------------
! Numerical Recipes subroutine.
!
! Contents:
! GAUSSJ
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Linear equation solution by Gauss-Jordan elimination, eq. (2.1.1).
! a(1:n,1:n) is an input matrix stored in an array of physical
! dimensions np by np. b(1:n,1:m) is an input matrix containing the m
! right-hand side vectors, stored in an array of physical dimensions np
! by mp. On output, a(1:n,1:n) is replaced by its matrix inverse, and
! b(1:n,1:m) is replaced by the corresponding set of solution vectors.
!-----------------------------------------------------------------------
	subroutine GAUSSJ(a,n,np,b,m,mp)

	implicit none
	integer m,mp,n,np,NMAX
	double precision a(np,np),b(np,mp)
! The largest anticipated value of n.
	parameter(NMAX=50)

	integer i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
	double precision big,dum,pivinv

	do j=1,n
	   ipiv(j)=0
	enddo
	do i=1,n
	   big=0d0
	   do j=1,n
	      if(ipiv(j).ne.1)then
		 do k=1,n
		    if(ipiv(k).eq.0)then
		       if(abs(a(j,k)).ge.big)then
			  big=abs(a(j,k))
			  irow=j
			  icol=k
		       endif
		    else if(ipiv(k).gt.1)then 
		       stop "GAUSSJ: Singular matrix! (1)"
		    endif
		 enddo
	      endif
	   enddo
	   ipiv(icol)=ipiv(icol)+1
	   if(irow.ne.icol)then
	      do l=1,n
		 dum=a(irow,l)
		 a(irow,l)=a(icol,l)
		 a(icol,l)=dum
	      enddo
	      do l=1,m
		 dum=b(irow,l)
		 b(irow,l)=b(icol,l)
		 b(icol,l)=dum
	      enddo
	   endif
	   indxr(i)=irow
	   indxc(i)=icol
	   if(a(icol,icol).eq.0d0) stop "GAUSSJ: Singular matrix! (2)"
	   pivinv=1d0/a(icol,icol)
	   a(icol,icol)=1d0
	   do l=1,n
	      a(icol,l)=a(icol,l)*pivinv
	   enddo
	   do l=1,m
	      b(icol,l)=b(icol,l)*pivinv
	   enddo
	   do ll=1,n
	      if(ll.ne.icol)then
		 dum=a(ll,icol)
		 a(ll,icol)=0d0
		 do l=1,n
		    a(ll,l)=a(ll,l)-a(icol,l)*dum
		 enddo
		 do l=1,m
		    b(ll,l)=b(ll,l)-b(icol,l)*dum
		 enddo
	      endif
	   enddo
	enddo

	do l=n,1,-1
	   if(indxr(l).ne.indxc(l))then
	      do k=1,n
		 dum=a(k,indxr(l))
		 a(k,indxr(l))=a(k,indxc(l))
		 a(k,indxc(l))=dum
	      enddo
	   endif
	enddo

	return
	end
