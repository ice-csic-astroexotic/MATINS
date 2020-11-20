!-----------------------------------------------------------------------
! Changes:
! wrt version 2:
! - Some simplifications.
! wrt version 1:
! - Substantial changes.
! - Printing of output moved to the main program.
! - A minimum size for dtb set.
! wrt version 0:
! - Code reformatted.
! - New timestep.h definition file created.
!
! To-do:
! - Improve the calculation of Courant velocities (cf. Daniele thesis
! pages 66-70). See CFL (Courant-Friedrich-Lux).
! - Compare with the original and check for inadvertent omissions.
!
! Contents:
! BTIMESTEP
!-----------------------------------------------------------------------
  subroutine BTIMESTEP(dtb,nits)
      use grid, only: fh, dtb_cour
      use input_params, only : kcour
      use grid, only: NANG, NP, KMAX, LMAX, rb, theta, bm
    
	implicit none
      real*8 dt, tbyear, tsnapad, tyear

! Arguments.
	integer nits
	real*8 dtb,j2(KMAX,LMAX)


! Definitions specific to the subroutine.
! Note that the save statement below is already implied by construction.
	integer, save :: icount=0
	integer i,j,k,l,maxv(2)
	real*8 dtc,cour(KMAX,LMAX),factor
	real*8 velmax,vwhist,vhall

!-----------------------------------------------------------------------
! Notes:
! Courant time is defined as dt = k*min(dr/vel).
! k is given in the input file (safe value for Hall is 1d-2).
! The maximum is evalued everywhere except at the outer boundary, where
! J is much higher.
! fh*|B| is in km^2/Myr, dr in km.
!-----------------------------------------------------------------------

! Calculation of dtb.
! dtbfix is read from the input file.
	dtc=0d0

! Initialize the array cour to some large value.
! The array cour is only used in this subroutine.
! Warning: The value of ia seems to be initialized by BEVOL, which comes
! after BTIMESTEP!
	cour=1d10
	   do l=2,LMAX
	      do k=2,KMAX
		 i=2*k-2
		 j=2*l-1
! Calculation of velocities.
! vhall is now calculated in a way that avoids division by zero.
		 vwhist = fh(j)*dsqrt(j2(k,l))
cc		 vhall=bm(i,j)*fh(j)**2
cc     &		 *(1d0/fh(j-1)-1d0/fh(j+1))/(rb(j+1)-rb(j-1))
		 vhall=bm(i,j)
     &		 *dabs((fh(j+1)-fh(j-1))/(rb(j+1)-rb(j-1)))

		 velmax=max(vwhist,vhall,
     &		 sqrt(vr_amb(i,j)**2+vth_amb(i,j)**2+vphi_amb(i,j)**2))

		 if(velmax.gt.0d0)then
		    cour(k,l)=min(rb(j+1)-rb(j-1),rb(j)*(theta(i+1)-theta(i-1)))
     &		    /(velmax*1.d-6)
		 endif
	      enddo
	   enddo

! Find the minimum value of the array cour and its location.
	   dtc=minval(cour)
	   maxv=minloc(cour)

	   factor=kcour

! Adjust the timestep.
! Warning: dtbmax was not defined anywhere and has been removed!
cc	   dtb=max(dtbfix,min(dtc*factor,dtbmax))
	   dtb=dtc*factor

! Check the value of dtb.
! Can be commented out to speed up execution.
	if(dtb.le.0d0)then
	   write(*,100)"BTIMESTEP: Error! Negative or zero dtb."
	   stop
	endif

! Set a minimum size for dtb wrt dt to avoid prolonged executions.
	if(dt/dtb.gt.1d4)then
	   dtb=dt/1d4
	   write(*,100)"BTIMESTEP: Warning! Courant time is very small, reset to: ",dtb
	endif

! Fixed number of iterations, but at least 2.
! dtb may need to be adjusted here to synchronize tbyear and tyear.
	nits=int(dt/dtb)+2

! Write out some diagnostic data.
	icount=icount+1
	if(icount.eq.1)then
	   open(2,file="diag/btimestep_count.txt")
	   write(*,100)"BTIMESTEP: File btimestep_count.txt written."
	   write(*,*)
	   write(2,100)"# Time step data."
	   write(2,200)"#","icount","nits","dtb","tbyear",
     &	   "dt","tyear","tsnapad","vwhist","vhall"
	else
	   open(2,file="diag/btimestep_count.txt",access="append")
	endif
	write(2,300)icount,nits,dtb,tbyear,dt,tyear,tsnapad,vwhist,vhall
	close(2)

! Format statement.
100	format(a)
200	format(a1,a9,a10,7a17)
300	format(2i10,7e17.8)

	return
	end
