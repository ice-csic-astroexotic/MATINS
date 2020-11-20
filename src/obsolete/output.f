      SUBROUTINE OUTPUT(icount,tyear,pc,sfluxb,tss,eta,omegatau,
     &		qcj,emnu,shearm,sigma_max,per,pdot,bindex,bpdip,poytot,
     &            joutot,joucortot)
      use constants
      
      implicit none
      include 'decl/dim2.h' 
      include 'decl/varsB.h'

      integer kslice,lslice
      parameter (kslice=kd/2, lslice=lmax-2)

      integer icount,i,j,k,l
      real*8 tyear
      real*8 per,pdot,bpdip,bindex,sd_age

      real*8 pc(kd,ld),  emnu(kd,ld)
      real*8 qcj(kd,ld), eta(kd,ld), omegatau(kd,ld)

      real*8 lumerg,lumph,teffi,tspole,tseq
      real*8 sfluxb(kd),tss(kd)
      real*8 phflux(kd)
      real*8 shearm(kd,ld), sigma_max(kd,ld)
      real*8 poytot,joutot,joucortot

! ============================
! 	MICROPHYSICS
! ============================

      if (tyear .eq. 0.d0) then
	open(unit=91,file='out/etaang.yg')
	open(unit=92,file='out/etarad.yg')
	open(unit=93,file='out/fh.yg')
        do l=1,lmax
	  write(93,*) rb(2*l-1),fh(2*l-1)
	enddo
	close(93)
        open(unit=99,file='out/omegatau.yg')

      else
	open(unit=91,file='out/etaang.yg',access='append')
	open(unit=92,file='out/etarad.yg',access='append')
        open(unit=99,file='out/omegatau.yg',access='append')
      endif

      write(91,80)'"Label= eta (',r(lslice),',th)'
      write(92,85)'"Label= eta (r,',z(kslice),')'
      write(99,85)'"Label=omtau(r,',z(kslice),')'

      write(91,*)'"Time=', tyear
      write(92,*)'"Time=', tyear
      write(99,*)'"Time=', tyear

      do k=1,kmax
	write(91,*) z(k),eta(k,lslice)
      enddo

      do l=2,lmax
        write(92,*) rb(2*l-1),eta(kslice,l)
	write(99,*) rb(2*l-1),omegatau(kslice,l)
      enddo

      write(91,*)
      write(92,*)
      write(99,*)

      close(91)
      close(92)
      close(99)

   80 format(a14,f5.2,a4)
   85 format(a15,f4.2,a1)
   88 format(a24,f4.2,a1)

c=====================================================================
c   calculate T_eff[K]
c   integrated luminosity, sfluxb in 10**40 erg/km**2/s, arl in km2
c=====================================================================


      phflux = PHFLUX_CONSTANT*UNIT_R**2*tss**3 ! photon flux in ph/km^2/s
      lumph = 0.d0
      lumerg  = 0.d0
      do k=2,kmax
        lumerg  = lumerg  + sfluxb(k)*
     &          2.d0*pi*r(lmax)**2*(dcos(z(k-1))-dcos(z(k)))
        lumph = lumph + phflux(k)*
     &          2.d0*pi*r(lmax)**2*(dcos(z(k-1))-dcos(z(k)))
      enddo
      lumerg=lumerg*1.d40
      teffi=(lumerg/(STEFAN_BOLTZMANN*4.d0*pi*(UNIT_R*r(lmax))**2))**(0.25d0)

      if (pdot .ne. 0.) sd_age=per/(pdot*T_YEAR)

! =======================================
! 	WRITING OUTPUTS
! =======================================

      open(unit=14,file='out/lastmaps.dat')
      if (tyear .eq. 0.d0) then
	open(unit=13,file='out/Tmap.dat')
	open(unit=35,file='out/Tp.yg')
	open(unit=36,file='out/Te.yg')
	open(unit=37,file='out/Tb.yg')
	open(unit=39,file='out/Ts.yg')
        open(unit=38,file='out/cool_curve.d')
	open(unit=45,file='out/qja.yg')
	open(unit=46,file='out/qjr.yg')
	open(unit=47,file='out/qnua.yg')
	open(unit=48,file='out/qnur.yg')
      else
	open(unit=13,file='out/Tmap.dat',access='append')
	open(unit=35,file='out/Tp.yg',access='append')
	open(unit=36,file='out/Te.yg',access='append')
	open(unit=37,file='out/Tb.yg',access='append')
	open(unit=39,file='out/Ts.yg',access='append')
        open(unit=38,file='out/cool_curve.d',access='append') 
	open(unit=45,file='out/qja.yg',access='append')
	open(unit=46,file='out/qjr.yg',access='append')
	open(unit=47,file='out/qnua.yg',access='append')
	open(unit=48,file='out/qnur.yg',access='append')
      endif

      write(38,210) tyear, lumerg, teffi, bpdip*UNIT_B, per,
     &   pdot, sd_age, pc(2,lmax)*UNIT_T, pc(kmax/2+1,lmax)*UNIT_T,
     &	 tss(2), tss(kmax/2+1), bindex, lumph
      close(38)

      write(*,220) 'COOLING',icount,'t=',tyear,'--Tp,Teq[K]=',
     &  pc(2,lmax)*UNIT_T,pc(kmax/2+1,lmax)*UNIT_T,'--Tsp,Tse[K]=',
     &  tss(2),tss(kmax/2+1),' Bpdip [G]=',bpdip*UNIT_B,' Lerg/ph=',lumerg,lumph

210   format(13(1x,1pe15.7))
220   format(a7,i7,a3,1pe10.3,2(a13,1p2e10.3),a6,1pe10.3,a9,2(1pe10.3))

      write(35,*)'"Time=', tyear
      write(35,*)'"Label= Tpole'
      write(36,*)'"Time=', tyear
      write(36,*)'"Label= Tequa'
      write(37,*)'"Time=', tyear
      write(37,*)'"Label= Tb'
      write(39,*)'"Time=', tyear
      write(39,*)'"Label= Tsur'
      write(45,*)'"Time=', tyear
      write(45,80)'"Label= qj  (',r(lslice),'th)'
      write(46,*)'"Time=', tyear
      write(46,85)'"Label= qj  (r,',theta(kslice),')'
      write(47,*)'"Time=', tyear
      write(47,80)'"Label= qnu (',r(lslice),'th)'
      write(48,*)'"Time=', tyear
      write(48,85)'"Label= qnu (r,',theta(kslice),')'

      write(14,*) tyear,per,pdot
      write(14,*) joutot,joucortot,poytot
      do i=0,nang+1
      do j=0,np+2
	write(14,*) br(i,j),bth(i,j),bphi(i,j),aphi(i,j)
      enddo
      enddo
      do k=1,kd
      do l=1,ld
	write(14,'(2f10.3,1pe10.3)') pc(k,l)
      enddo
      enddo
      close(14)

      if (bpdip .ge. 1.d-2) then
      do k=2,kmax
      do l=2,lmax
	write(13,'(2f10.3,3(1pe10.3))') theta(2*k-2),rb(2*l-1),1.d8*pc(k,l),
     &		shearm(k,l),sigma_max(k,l)
      enddo
      enddo
      endif
      close(13)

      do l=2,lmax
        write(35,*) rb(2*l-1),1.d8*pc(2,l)
        write(36,*) rb(2*l-1),1.d8*pc(kmax/2,l)
        write(46,*) rb(2*l-1),qcj(kslice,l)
        write(48,*) rb(2*l-1),emnu(kslice,l)
      enddo

      do k=2,kmax
        write(37,*) theta(2*k-2),dlog10(UNIT_T*pc(k,lmax))
        write(39,*) theta(2*k-2),dlog10(tss(k))
        write(45,*) theta(2*k-2),qcj(k,lslice)
        write(47,*) theta(2*k-2),emnu(k,lslice)
      enddo

      write(35,*)
      write(36,*)
      write(37,*)
      write(39,*)
      write(45,*)
      write(46,*)
      write(47,*)
      write(48,*)

      close(35)
      close(36)
      close(37)
      close(39)
      close(45)
      close(46)
      close(47)
      close(48)
 
      return
      end
