!-----------------------------------------------------------------------
! Test array communication to subprograms in Fortran.
!
! Fortran stores arrays by column: A11,A21,A31,...,A12,A22,A32,...
!-----------------------------------------------------------------------
	implicit none
	integer i,j,ni,nj
	parameter(ni=5,nj=5)
	integer A(ni,nj)

	do i=1,ni
	   do j=1,nj
	      A(i,j)=(i-1)*nj+j
	   enddo
	enddo

	do i=1,ni
	   write(*,*)(A(i,j),j=1,nj)
	enddo

	call ARRAY_CHECK(ni-2,nj-2,A)

	call ARRAY_CHECK(ni-2,nj-2,A(2:ni-1,2:nj-1))

	call ARRAY_CHECK(ni-2,nj-2,A(3:ni,3:nj))

	stop
	end
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
	subroutine ARRAY_CHECK(ni,nj,A)

	implicit none
	integer i,j,ni,nj
	integer A(ni,nj)

	write(*,*)
	do i=1,ni
	   write(*,*)(A(i,j),j=1,nj)
	enddo

	return
	end
