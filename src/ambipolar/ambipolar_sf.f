      subroutine ambipolar_sf(lcore, t8, lambda, lambda_DU, tau_pn)
c       -------------------------------
      use structure, only: rho,yn,yp,tccru,tcn,tcp
      implicit none
!-----------------------------------------------------------------------
! Dimensions of the grid.
!-----------------------------------------------------------------------
      integer m, j, lcore
      include '../decl/dim2.h'
      real*8, dimension(lmax) :: lambda,lambda_DU,tau_pn
C
      real*8 rho_nuc, mp, nbaryon, t8, pi
      real*8 rho0, rhomed, ynmed, ypmed, tcnmed,tcpmed

      real*8 kFn, kFp, vg_p, vg_n, taun, taup, v1, v2
      real*8 Rp_1, Rp_2, RD_sf, RM_sf
      real*8 RRp_1p, RRp_2p, RRD_A, RRMp_pA, RRMn_pA
      real*8 RRp_1n, RRp_2n, RRD_B, RRMn_nB, RRMp_nB
      real*8 RRp_1, RRp_2, RRD_BA, RRMn_BA, RRMp_BA
      real*8 Kp_1,Kp_2,Sp_1,Sp_2, R_pn_0, R_pn, Rsf

      rho_nuc=2.8d14
      mp = 1.672d-24  ! proton mass
      pi = dacos(-1.d0)

      do m=1,lcore
        rho0=rho(2*m)/rho_nuc
        lambda(m) = 5d27*t8**6*rho0**(2d0/3d0)
        lambda_DU(m)=3.5d36*t8**4*rho0**(1d0/3d0)
        tau_pn(m)=rho0**(1d0/3d0)/(4.7d16*t8**2)
cold        beta(m)=1.02d34*rho0**(4d0/3d0)/t8**2
      enddo
C-----------------------------------------------------------
C    SUPERFLUID/SUPERCON CORRECTIONS
C-----------------------------------------------------------

      do m=1,lcore
        j = 2*m
        tcnmed=tcn(j)
	      tcpmed=tcp(j)
        taun=t8/dmax1(tcnmed,0.9d0*t8)
        taup=t8/dmax1(tcpmed,0.9d0*t8)

        if (taup.ge.1d0.and.taun.ge.1d0) then ! p norm, n norm 
          Rp_1 = 1d0
          Rp_2 = 1d0
          RD_sf = 1d0
          RM_sf = 1d0
        elseif (taup.lt.1d0.and.taun.ge.1d0) then ! p sup,  n normal
          vg_p=v1(taup)   !1S_0 (nB)  
          Rp_1  = RRp_1p(vg_p)
          Rp_2  = RRp_2p(vg_p)
          RD_sf = RRD_A(vg_p)
          RM_sf = max(RRMp_pA(vg_p),RRMn_pA(vg_p))
        elseif (taup.ge.1d0.and.taun.lt.1d0) then ! p norm, n sup
          vg_n=v2(taun)   !3P_2 (nB)  
          Rp_1  = RRp_1n(vg_n)
          Rp_2  = RRp_2n(vg_n)
          RD_sf = RRD_B(vg_n)
          RM_sf = max(RRMn_nB(vg_n),RRMp_nB(vg_n))
        elseif (taup.lt.1d0.and.taun.lt.1d0) then ! p and n sup
          vg_p=v1(taup)   !1S_0 (nB)  
          vg_n=v2(taun)   !3P_2 (nB)  
          Rp_1  = RRp_1(vg_n,vg_p)
          Rp_2  = RRp_2(vg_n,vg_p)
          RD_sf = RRD_BA(vg_n,vg_p)
          RM_sf = max(RRMn_BA(vg_n,vg_p),RRMp_BA(vg_n,vg_p))
        endif
        rhomed=rho(j)
        ynmed=yn(j)
	      ypmed=yp(j)
        nbaryon = 1.d-39*rhomed/mp  !! in fm-3
        kFn = (3d0*pi**2*nbaryon*ynmed)**(1.d0/3.d0) !! in fm-1
        kFp = (3d0*pi**2*nbaryon*ypmed)**(1.d0/3.d0) !! in fm-1

        R_pn_0 = Sp_2(kFn,kFp)  * Kp_2(kFn,kFp)
     &            + Kp_1(kFn,kFp)  * Sp_1(kFn,kFp)

        R_pn = Sp_2(kFn,kFp)  * Kp_2(kFn,kFp) * Rp_2
     &             + 0.5d0 * Kp_1(kFn,kFp) * Sp_1(kFn,kFp)
     &             * ( 3d0 * Rp_1 - Rp_2 )
        R_pn = R_pn/R_pn_0

        lambda(m) = lambda(m)*RD_sf
        lambda_DU(m) = lambda_DU(m)*RM_sf
        tau_pn(m)=tau_pn(m)/R_pn
      enddo
C-----------------------------------------------------------
C-----------------------------------------------------------

      return
      end


