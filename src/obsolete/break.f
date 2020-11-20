C------------------------------------------------------------
      SUBROUTINE BREAKING(dtb,time,shearm,smax,iwrite)
C------------------------------------------------------------
      implicit none
      include '../decl/dim2.h'
      include '../decl/grid.h'
      include '../decl/varsB.h'
      integer kbr,lbr,k,l,i,j,iwrite

      real*8 fac, sumen, sumvol, eps, time, tign, dtb
      real*8 shearm(kd,ld), smax(kd,ld)
      real*8 fr0(kd,lcore:lmax), fth0(kd,lcore:lmax),
     &		fphi0(kd,lcore:lmax)
      real*8 freq(kd,ld), ftheq(kd,ld), fphieq(kd,ld)
      real*8 str(kd,lcore:lmax), stth(kd,lcore:lmax),
     &		stphi(kd,lcore:lmax)

      integer indx(kd,lcore:lmax)
      logical breakr,breakth,breakphi, breaktot

      save freq,ftheq,fphieq
      parameter(tign=0.d0)

C_-------------------------------------------------------------------------
C_-------------------------------------------------------------------------

      eps = 0.9d0
      fac = 1.d24/(4.d0*PI)   ! B12 --> Gauss
      do k=2,kmax
      do l=lcore+1,lmax
        i=2*k-2
        j=2*l-1

        fr0(k,l)   = fac*bth(i,j)*bphi(i,j)
        fth0(k,l)  = fac*bphi(i,j)*br(i,j)
        fphi0(k,l) = fac*br(i,j)*bth(i,j)
	if (time.le.tign) then
	  freq(k,l)   = fr0(k,l)    
	  ftheq(k,l)  = fth0(k,l)    
	  fphieq(k,l) = fphi0(k,l)
	endif
      enddo
      enddo

      if (time.eq.0.d0) then
	open(unit=50,file='outb/break.d',status='replace')
	open(unit=51,file='outb/mapbreak.d',status='replace')
	close(50)
	close(51)
      endif

      str  =0.d0
      stth =0.d0
      stphi=0.d0

      do k=2,kmax
      do l=lcore+1,lmax
        str(k,l)  =    dabs(fr0(k,l)-freq(k,l))/smax(k,l)
        stth(k,l) =  dabs(fth0(k,l)-ftheq(k,l))/smax(k,l)
        stphi(k,l)=dabs(fphi0(k,l)-fphieq(k,l))/smax(k,l)
      enddo
      enddo

      if (iwrite .eq. 1 .and. time .gt. tign) then
	open(unit=51,file='outb/mapbreak.d',access='append')
	write(51,'(1e12.3)') time
	do k=2,kmax
	do l=lcore+1,lmax
	  write(51,'(5f10.3)') zc(k),rc(l),str(k,l),stth(k,l),stphi(k,l)

	if ((str(k,l) .ne. str(k,l)) .or.
     & (stth(k,l) .ne. stth(k,l)) .or.
     &	(stphi(k,l) .ne. stphi(k,l))) then
	  print*,'FAIL break.f'
	  print*,k,l,smax(k,l),str(k,l),stth(k,l),stphi(k,l)
	  print*,k,l,dabs(freq(k,l)),time
	  stop
	endif

	enddo
	enddo
	close(51)
      endif

      do k=2,kmax
      do l=lcore+1,lmax
        breakr   = (str(k,l)   .gt. eps)
        breakth  = (stth(k,l)  .gt. eps)
        breakphi = (stphi(k,l) .gt. eps)
        if (breakr.or.breakth.or.breakphi) then
          indx(k,l)=1
        else
          indx(k,l)=0
        endif
      enddo
      enddo

      breaktot=.false.
      do k=2,kmax
      do l=lcore+1,lmax
        breakr   = (str(k,l)   .gt. 1.d0)
        breakth  = (stth(k,l)  .gt. 1.d0)
        breakphi = (stphi(k,l) .gt. 1.d0)

        if (breakr.or.breakth.or.breakphi) then
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
	do l=lcore+1,lmax
	  if(indx(k,l).eq.1) then
	    sumen = sumen + ((str(k,l)*smax(k,l))**2 +
     &	      (stth(k,l)*smax(k,l))**2+(stphi(k,l)*smax(k,l))**2)
     &        /shearm(k,l)*(vol(2*k-2,2*l-1)*1.d15)
	    sumvol=sumvol+(vol(2*k-2,2*l-1)*1.d15)
	    freq(k,l)=fr0(k,l)
	    ftheq(k,l)=fth0(k,l)
	    fphieq(k,l)=fphi0(k,l)
	  endif  
	enddo
	enddo

	if (time.gt.0.) then
	  open(unit=50,file='outb/break.d',access='append')
	  write(*,444) 'CRUST BREAKS, ENERGY=',
     &       sumen,time,zc(kbr),rc(lbr),sumvol, dtb

	  write(50,445) time,sumen/1.d40,zc(kbr),rc(lbr),sumvol,dtb
	  close(50)
	endif
      endif

  444 format(A24, 2x, es11.3, es14.6, 2f8.3, 3x, e10.3, 7e12.3) 
  445 format(es12.6, 3x, es9.3, 2f8.3, 3x, 2e10.3) 

      return
      end
