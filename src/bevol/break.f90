!> In this subroutine, we aim to calculate 
!> The crustal failures, under the approach of Perna & Pons 2011
!>
!> @authors:
!> Jose Pons
!> Clara Dehman


module magnetic_stress

  use utils, only: get_free_unit
  use constants, only: UNIT_B, UNIT_R, PI
  use grid, only: br, bth, bphi, KMAX, LMAX, benu, vol
  use grid, only: shearModulus, shearMaximum, jcore
  use grid, only: freq, ftheq, fphieq
  use grid, only: theta, rb
  
  contains

    subroutine breaking(dtb,time)
    implicit none
    integer kbr, lbr, k, l, i, j
    integer :: unit_1, unit_2 
    real*8, parameter :: tign=0.d0, eps = 0.9d0   ! tign is the time to be ignored
    real*8, parameter :: nwrite=1.d0
    real*8 sumen, sumvol
    real*8 fac
    real*8, intent(in) :: time, dtb
    real*8, dimension(kmax,jcore/2:lmax) :: fr0, fth0, fphi0
    real*8, dimension(kmax,jcore/2:lmax) :: str, stth, stphi
    integer indx(kmax,jcore/2:lmax)
    logical breakr, breakth, breakphi, breaktot
    character(len=68), parameter :: break_format = &
     &   "(es12.6, 3x, es9.3, 2f8.3, 3x, 2e10.3)"
      character(len=68), parameter :: CRUST_BREAKS_ENERGY = &
     &   "(A24, 2x, es11.3, es14.6, 2f8.3, 3x, e10.3, 7e12.3)"
      character(len=16), parameter :: time_format = &
     &   "(1e12.3)"
      character(len=16), parameter :: mapbreak_format = &
     &   "(5f10.3)"


!>----------------- output files -----------------------------------------
    unit_1 = get_free_unit()
    unit_2 = get_free_unit()


    fac = UNIT_B**2.d0/(4.d0*PI)   ! B12 --> Gauss
    do k=2,kmax
      do l=jcore/2+1,lmax
        i=2*k-2
        j=2*l-1
        fr0(k,l)   = fac*bth(i,j)*bphi(i,j)
        fth0(k,l)  = fac*bphi(i,j)*br(i,j)
        fphi0(k,l) = fac*br(i,j)*bth(i,j)
        if (time <= tign) then
          freq(k,l)   = fr0(k,l)    
          ftheq(k,l)  = fth0(k,l)    
          fphieq(k,l) = fphi0(k,l)
        endif
      end do
    end do


    if (time == 0.d0) then
      open(unit=unit_1,file='outb/break.d',status='replace')
      open(unit=unit_2,file='outb/mapbreak.d',status='replace')
      close(unit_1)
      close(unit_2)
    endif

    str  =0.d0
    stth =0.d0
    stphi=0.d0
    do k=2,kmax
      do l=jcore/2+1,lmax
        str(k,l)  =    dabs(fr0(k,l)-freq(k,l))/shearMaximum(k,l)
        stth(k,l) =  dabs(fth0(k,l)-ftheq(k,l))/shearMaximum(k,l)
        stphi(k,l)=dabs(fphi0(k,l)-fphieq(k,l))/shearMaximum(k,l)
      enddo
    enddo

    if (modulo(time,nwrite) == 0.d0 .and. time > tign) then
      open(unit=unit_2,file='outb/mapbreak.d',access='append')
      write(unit_2,time_format) time
      do k=2,kmax
        do l=jcore/2,lmax
          write(unit_2,mapbreak_format) theta(2*k-2),rb(2*l-1),str(k,l),stth(k,l),stphi(k,l)
        enddo
      enddo
      close(unit_2)
    endif


    do k=2,kmax
      do l=jcore/2+1,lmax
        breakr   = (str(k,l) > eps)
        breakth  = (stth(k,l) > eps)
        breakphi = (stphi(k,l) > eps)
        if (breakr.or.breakth.or.breakphi) then
          indx(k,l)=1
        else
          indx(k,l)=0
        endif
      enddo
    enddo

    breaktot=.false.
    do k=2,kmax
      do l=jcore/2+1,lmax
        breakr   = (str(k,l) > 1.d0)
        breakth  = (stth(k,l) > 1.d0)
        breakphi = (stphi(k,l) > 1.d0)
        if (breakr .or. breakth .or. breakphi) then
          breaktot=.true.
          kbr=k
          lbr=l
          goto 40
        endif
      enddo
    enddo

 40 continue
    if (breaktot) then
      sumen = 0.d0
      sumvol = 0.d0
      do k=2,kmax
        do l=jcore/2+1,lmax
          if(indx(k,l) == 1) then
            sumen = sumen + ((min(str(k,l),1.d0)*shearMaximum(k,l))**2+(min(stth(k,l),1.d0)*shearMaximum(k,l))**2+ &
             &    (min(stphi(k,l),1.d0)*shearMaximum(k,l))**2)/shearModulus(k,l)*(vol(2*k-2,2*l-1)*UNIT_R**3.d0)
            sumvol=sumvol+(vol(2*k-2,2*l-1)*UNIT_R**3.d0)
            freq(k,l)=fr0(k,l)
            ftheq(k,l)=fth0(k,l)
            fphieq(k,l)=fphi0(k,l)
          endif  
        enddo
      enddo
  
      if (time > 0.d0) then
        open(unit=unit_1,file='outb/break.d',access='append')
        write(unit_1,break_format) time,sumen,theta(2*kbr-2),rb(2*lbr-1),sumvol,dtb
        close(unit_1)
      endif
    endif

    end subroutine breaking

end module magnetic_stress
