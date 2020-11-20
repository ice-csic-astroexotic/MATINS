       subroutine cphonon_Chu(t,zh,ah,ye,rho,xh,yn,Cion,k_tot,L_ie)
c       Chugunov & Haensel 2007
       implicit none
       double precision t,rho,nbaryon
       double precision ye,ah,zh,xh,yn
       double precision gamma,nions
       double precision xi
       double precision k_ii,k_ie,k_tot,kb,nel,k_F,k_Fless
       double precision clight,mu,me,hbar,pi
       double precision cs,q_BZ,Cion,L_ie,v_F,v_F2,eexp,u1,u2
       double precision w_DW,yy,lambda_ie,XEE1,XEE1_1,expintt
       double precision arg1, arg2,FF
       double precision xe
       double precision rho6,k_0,k_ast,beta,wp,Tp,theta,ai,e13
       double precision aprime, a2prime

       external FF, expintt
       kb=1.381d-16 !erg/K
       me=9.11d-28 !gr
       mu=1.66d-24 !gr
       clight=2.9979d10 !cm/s
       hbar=1.0546d-27  !erg*s
       pi=3.141592653589793d0
       rho6= rho/1.0d6 
       e13=1.d0/3.d0

       xe=1.009d0*((rho6*ye)**(e13)) !electron relat. parameter

       nbaryon=rho/1.66d-24   !1/cm3
       nions=nbaryon*xh/ah    !1/cm3
       gamma = 2.275d5*zh**2*(rho*xh/ah)**(e13)/t

      IF (gamma.le.1.d0) THEN  ! a gas
        k_ii=0.d0
        k_ie=0.d0
        k_tot=0.d0
        return
      ELSE    !a solid or liquid
c      alat=(ah*mu/rho)**(1./3.)
c     ---------phonon-phonon Kii ---------------------------
       k_ast= 0.4d0
       beta = e13
       ai= (4.d0*pi*nions/3.d0)**(-e13) !ion sphere radius (cm)
       a2prime=yn*ah/xh    !check with a2prime=yn*ah
       aprime=ah + a2prime
       Tp = 7.832d6*dsqrt(rho6*zh**2/aprime/ah) !plasma T in K
       theta = Tp/t 
       wp = kb*Tp/hbar     !plasma freq. in 1/s
       k_0 = kb*wp*nions*ai**2
       xi=0.7d0    
       k_ii=k_0*
     - dsqrt(k_ast**2+(xi**3*gamma/2.85d0)**2*dexp(2.d0*beta*theta))
c     ---------phonon-electron Kie----------------------------------
       q_BZ=(6*pi**2*nions)**e13  !q of Brioulin zone in 1/cm
c       cs=wp/3.d0/q_BZ            !sound speed-group velocity
       cs=0.7d0*wp/q_BZ  !Itoh
       
       nel =nbaryon*ye  !1/cm3
       k_F = (3*pi**2*nel)**(e13) !1/cm
       k_Fless =2.d0*k_F/q_BZ          !dimensionless y^{-1}
       v_F=clight*xe/dsqrt(1.d0+xe**2) !electron Fermi veloc. 
       v_F2=(v_F/clight)**2            !dimensionless
       u1=2.8d0
       u2=13.d0
       w_DW=1.683d0*dsqrt(xe/ah/zh)
     -  *(0.5d0*u1*dexp(-9.1d0/theta)+u2/theta)
       yy=1.d0/k_Fless 
       arg1=w_DW*yy**2
       arg2=w_DW      
       eexp=(dexp(-arg1)-dexp(-arg2))/arg2
       XEE1=expintt(1,arg1)
       XEE1_1=expintt(1,arg2)
       lambda_ie = 0.5d0*(XEE1-XEE1_1-v_F2*eexp)
       L_ie=(320.d0*ai/(1.d0+xe**2)/lambda_ie)*(26.d0/zh)
     -  *(FF(theta)/0.01d0)*dsqrt(ah*rho6/aprime)   !mean free path  eq.43     
c       Cion=ccvfion(t,rho,xh,ah,zh)
c      
       k_ie=e13*Cion*cs*L_ie    !*1.d3  if  we put a factor 1e3, their results are 
c       similar to our phonons.  
c      -------------Total-------------------------------
       k_tot=1.d0/(1.d0/k_ii+1.d0/k_ie)

      ENDIF

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccc      
      double precision FUNCTION FF(theta)
      implicit none
      double precision theta

      FF=0.014d0+0.03d0/(dexp(theta/5.d0)+1.d0)
      return
      end
cccccccccccccccccccccccccccccccccc
      double precision FUNCTION expintt(n,xx)
      INTEGER n,MAXIT
      double precision x,EPS,FPMIN,EULER
      PARAMETER (MAXIT=100,EPS=1.e-7,FPMIN=1.e-30,EULER=.5772156649)
      INTEGER i,ii,nm1
      double precision a,b,c,d,del,fact,h,psi,xx
      nm1=n-1
      x=xx
c      print*, n, x
      if(n.lt.0.or.x.lt.0..or.(x.eq.0..and.(n.eq.0.or.n.eq.1)))then
        write(*,*) "<warning>[CPHONON_CHU] Bad arguments in expint"
      else if(n.eq.0)then
        expintt=exp(-x)/x
      else if(x.eq.0.)then
        expintt=1./nm1
      else if(x.gt.1.)then
        b=x+n
        c=1./FPMIN
        d=1./b
        h=d
        do 11 i=1,MAXIT
          a=-i*(nm1+i)
          b=b+2.
          d=1./(a*d+b)
          c=b+a/c
          del=c*d
          h=h*del
          if(abs(del-1.).lt.EPS)then
            expintt=h*exp(-x)
            return
          endif
11      continue
        write(*,*) "<warning>[CPHONON_CHU] continued fraction failed in expint"
      else
        if(nm1.ne.0)then
          expintt=1./nm1
        else
         expintt=-log(x)-EULER
        endif
        fact=1.
        do 13 i=1,MAXIT
          fact=-fact*x/i
          if(i.ne.nm1)then
            del=-fact/(i-nm1)
          else
            psi=-EULER
            do 12 ii=1,nm1
              psi=psi+1./ii
12          continue
            del=fact*(-log(x)+psi)
          endif
          expintt=expintt+del
          if(abs(del).lt.abs(expintt)*EPS) return
13      continue
        write(*,*) "<warning>[CPHONON_CHU] Series failed in expint"
      endif
      return
      END

      

