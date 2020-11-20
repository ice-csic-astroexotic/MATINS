c  this subroutine calculates the radiative conductivities 
c  N.A. Silant'ev and  D.G. Yakovlev Ap&SS 71(1980) 45-50
c  for the bound-free process M.E. Schaaf A&A 227(1990) 61-70

       subroutine ykurrad(t,ye,rhoin,Bp,condphp,condpht)
       implicit none
       real*8 t,rhoin,condph,condpht,condphp,ne
       real*8 radkap,Bpnorm,ye,sigma,b,z1,z2
       real*8 radkapff,radkapes,f
       real*8 thompson
       real*8 condphof,condphoe,ap,at
       real*8 Bp
       
       ne=rhoin*ye/1.66d-24  !cgs
       sigma=5.67d-5  !cgs
       thompson=6.650643d-25 !cgs
       Bpnorm=Bp*1.d12/4.414d13
c
c free-free radiative opacity
       condphof=0.72d0*(t/1.d6)**(6.5d0)/(rhoin/1.d6)**2
       radkapff=16.d0*sigma*t**3/3.d0/rhoin/condphof  !cgs
c e-scattering radiative opacity
       radkapes=ne*thompson/rhoin
       condphoe=16.d0*sigma*t**3/3.d0/ne/thompson
c opacity for the two process
       radkap=radkapff+radkapes

c contribution of two process
       f=condphoe/(condphoe+condphof)
c total radiative conductivity
       condph=condphoe*condphof/(condphoe+condphof)
      
c parallel and opacity      

       b = (5.93014d+9*Bpnorm)/t
       call yakac2(b,f,z1,z2)
       ap = 1.d0/z1
       at = 1.d0/z2
       condphp=condph/ap
       condpht=condph/at

       return
       end

C--------------------------------------------------------------
      subroutine yakac2(bin,fin,zout1,zout2)
C--------------------------------------------------------------
c
      real z1(11,27),z2(11,27),f(11),b(27),sfin,binlog
      real*8 bin,db,df,fin,omdb,omdf,zlog1,zout1,zlog2,zout2

      integer ib, if, i, ibp, ifp
      save z1,z2,f,b
c
      zout1 = 1.
      zout2 = 1.
      if(bin.lt.1.) return
c
      if(fin.ge.0.and.fin.le.1.) go to 20
      write(66,21) fin
      write(6 ,21) fin
   21 format(' illeagal value of f, ',1pg11.4,', in yakac')
   20 continue
      sfin = sngl(fin)
      do 22 i = 1,10
      if = i
      ifp = i + 1
      if(f(if).le.sfin.and.sfin.lt.f(ifp)) go to 23
   22 continue
   23 continue
      df = (sfin - f(if))/(f(ifp) - f(if))
      omdf = 1.d0 - df
c
      binlog = sngl(dlog10(bin))
      ib = 26
      ibp = 27
      if(binlog.ge.b(27)) go to 26
      do 24 i = 1,26
      ib = i
      ibp = ib + 1
      if(b(ib).le.binlog.and.binlog.lt.b(ibp)) go to 25
   24 continue
   25 continue
   26 continue
c
      db = (binlog - b(ib))/(b(ibp) - b(ib))
      omdb = 1.d0 - db

      zlog1 =
     $       df*(db*z1(ifp,ibp) + omdb*z1(ifp,ib)) +
     $     omdf*(db*z1(if,ibp) + omdb*z1(if,ib))
      zout1 = 10.d0**zlog1

      zlog2 =
     $       df*(db*z2(ifp,ibp) + omdb*z2(ifp,ib)) +
     $     omdf*(db*z2(if,ibp) + omdb*z2(if,ib))
      zout2 = 10.d0**zlog2

      return
c
      entry yakitl
c
      open(23,file='yakdat2')   !,status='readonly')
c
      do 1 ib = 1,27
      read(23,*) b(ib)
      if(b(ib).gt.1.d3) pause
    1 continue

      do 2 if = 1,11
      read(23,*) f(if)
c      write(*,*)'f(if)',if,f(if)
      do 4 ib = 1,27
      read(23,*) z1(if,ib),z2(if,ib)
    4 continue
    2 continue
c
c            store b and z as logs (base 10)
c            f stored as read (no log)
      do 56 ib = 1,27
      b(ib) = alog10(b(ib))
      do 57 if = 1,11
      z1(if,ib) = alog10(z1(if,ib))
      z2(if,ib) = alog10(z2(if,ib))
   57 continue
   56 continue
c
      close(23)
c
      return
c  finish yakac2
      end
