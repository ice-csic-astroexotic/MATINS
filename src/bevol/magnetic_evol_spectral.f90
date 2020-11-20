! This subroutine calculates the structure of magnetic field with spectral methods
!---------------------------------------------------------------
!      phi -   poloidal field 
!      psi -   toroidal field 
!---------------------------------------------------------------  
subroutine magnetic_evol_spectral(dtb_myr)

  ! module imports -------------------------------------------------------------
  use grid, only : rb, cth, etab, enable_hall_effect
  use grid, only : np, nang, ng, nleg, belam, jevol
  use constants, only: PI
  use legpol, only: plnx, xg, wg, frel, fun_interp
  use initial_magnetic, only: potentials_to_b
  use initial_magnetic, only: phi, psi, phiold, psiold

  implicit none

  ! subroutine arguments -------------------------------------------------------
  real*8, intent(in) :: dtb_myr

  integer, parameter :: leg_filter = 20
  ! local variables ------------------------------------------------------------
  integer :: i, j, n
  real*8, dimension(0:np+2,nleg) :: phip, psip, phi1
  real*8, dimension(nang) :: eta_slice
  real*8, dimension(0:np+2,0:nleg) :: eta_spe
  real*8, dimension(0:np+2) :: drb
  real*8, dimension(nleg,nleg,np) :: aa, bb, cc
  real*8, dimension(nleg,np) :: rr, uu
  real*8, dimension(np,nleg) :: spol, stor
  real*8 :: fun, sum

  ! ----------------------------------------------------------------------------

  phiold = phi
  psiold = psi
  phip = 0d0
  psip = 0d0
  
  ! Definition of drb (similar to lr, but with no relativistic factors included)
  drb(0) = 0.5d0*rb(1)
  drb(np+2) = 0.5d0*(rb(np+2) - rb(np+1))
  do j = 1,np+1
    drb(j) = 0.5d0*(rb(j+1) - rb(j-1))

    ! Spectral decomposition of the magnetic diffusivity for each radial layer
    ! eta = sum eta_l P_l(mu)
    ! eta_l = (l + 1/2)*int(P_l eta(mu) dmu)

    eta_slice = etab(1:nang,j) 
    do n = 0,nleg
      sum = 0d0
      do i = 1,ng
        call fun_interp(nang,cth(1:nang),eta_slice,xg(i),fun)
        sum = sum + wg(i) * plnx(i,n) * fun
      end do
      ! Reset to zero very small values of the integral
      if (dabs(sum) < 1d-10) sum = 0d0
      eta_spe(j,n) = ( n + 0.5d0 )*sum
    end do
  end do

  ! calculating derivatives
  do j=jevol,np
    phip(j,:) = (phi(j+1,:)-phi(j-1,:))/(2d0*drb(j))
    psip(j,:) = (psi(j+1,:)-psi(j-1,:))/(2d0*drb(j))
    if (enable_hall_effect .eqv. .true.) then
      do n=1,nleg
      phi1(j,n) = (phi(j+1,n)-2.d0*phi(j,n)+phi(j-1,n))/(belam(j)*drb(j))**2  &
     &     -n*(n+1)*phi(j,n)/rb(j)**2
       end do
      endif
  end do

  ! non-linear terms (Hall)
  if (enable_hall_effect .eqv. .true.) then

    phi1(0,:) = 0.d0
    ! phi1(0,:) = -shat(0)*dnm(0,:)
    phi1(np+1:np+2,:) = 0.d0
 
    call hallterm(phip,psip,phi1,spol,stor)
  else
    spol=0d0 
    stor=0d0
  endif

  ! implicit time advance poloidal field
  call pol_matrix_inversion(dtb_myr,drb,eta_spe,aa,bb,cc)

  do n=1,nleg
    do j=1,np
      rr(n,j) = phiold(j,n) + dtb_myr*spol(j,n)
    end do
  end do

  call triblocb(np,nleg,aa,bb,cc,rr,uu)
  do n=1,nleg
    do j=1,np
      phi(j,n) = uu(n,j)
    end do
  end do

  ! boundary conditions for Phi
  phi(0,:) = 0.d0
  phip(0,:) = (phi(1,:)-phi(0,:))/drb(0)

  do n = 1,nleg
    phip(np+1,n) = - frel(n)*n*phi(np,n)/rb(np)
    phi(np+1,n) = phi(np,n) + drb(np+1)*phip(np+1,n)
    phip(np+2,n) = - frel(n)*n*phi(np+1,n)/rb(np+1)
    phi(np+2,n) = phi(np+1,n) + drb(np+2)*phip(np+2,n)
  end do

  ! implicit time advance toroidal field
  call tor_matrix_inversion(dtb_myr,drb,eta_spe,aa,bb,cc)

  do n=1,nleg
    do j=1,np
      rr(n,j) = psiold(j,n) + dtb_myr*stor(j,n)
    end do
  end do

  call triblocb(np,nleg,aa,bb,cc,rr,uu)

  do n=1,nleg
    do j=1,np
      psi(j,n) = uu(n,j)
    end do
  end do
  
  ! Boundary conditions for toroidal field
  psi(np:np+2,:) = 0.d0
  psip(np:np+2,:) = 0d0
  psi(np-1,:)=0.5d0*psi(np-2,:)

  ! filtering and cleaning
  do n=leg_filter,nleg
    phi(:,n) = phi(:,n)*dexp(-dble(n)/dble(leg_filter))
    psi(:,n) = psi(:,n)*dexp(-dble(n)/dble(leg_filter))
  end do
  do n=1,nleg
    do j=1,np
      if (dabs(phi(j,n)) < 1.d-30) phi(j,n) = 0.d0
      if (dabs(psi(j,n)) < 1.d-30) psi(j,n) = 0.d0
    end do
  end do

  call potentials_to_b
 
end subroutine magnetic_evol_spectral


subroutine hallterm(phip,psip,phi1,s1,s2)

  use initial_magnetic, only: phi, psi, bcg
  use grid, only: fh, np, rb
  use legpol, only: nleg
  use constants, only: PI
  implicit none

  real*8, dimension(0:np+2,nleg), intent(in) :: phip, psip, phi1
  real*8, dimension(np,nleg), intent(out) :: s1, s2

  integer :: i, j, n, k, kp, gg

  real*8, dimension(0:np+2,nleg) :: etheta, dnm, cnm, cnm2
  real*8 :: fac, coef2, coef3, i2, i3, bcf, icf

  ! Initialize values
  etheta = 0.d0
  dnm = 0.d0
  cnm = 0.d0
  cnm2 = 0.d0

  ! calculation of non-linear terms (Hall)
  do n=1,nleg
    do k=1,nleg
      do j=max(0,(n-nleg+k+1)/2),min(k,n)
        kp=k+n-2*j
        if (kp.gt.0) then
          gg = (kp+k+n)/2
          bcf = bcg(gg-kp)*bcg(gg-k)*bcg(gg-n)/bcg(gg)
          icf=dsqrt((2*k+1)*(2*kp+1)*(2*n+1)/(4*PI))/(2*gg+1)*bcf
          coef2 = 0.5d0*kp*(kp+1)* &
     &          (n*(n+1)-kp*(kp+1)+k*(k+1))/dble(n*(n+1))
          coef3 = 0.5d0*(-n*(n+1)+kp*(kp+1)+k*(k+1))

          i2 = coef2*icf
          i3 = coef3*icf

          ! i2 and i3 could be stored, calculated only once.

          do i=0,np+2
            dnm(i,n) = dnm(i,n)  &
     &          + i2*(phip(i,k)*psi(i,kp) - psip(i,k)*phi(i,kp))
            etheta(i,n)  = etheta(i,n)  &  
     &          + i2*(psi(i,k)*psi(i,kp)+phi1(i,k)*phi(i,kp))
            cnm(i,n) = cnm(i,n)   &
     &          + i3*(psi(i,k)*psip(i,kp) + phip(i,k)*phi1(i,kp))
            cnm2(i,n) = cnm2(i,n)  &
     &          + i3*(psi(i,kp)*psip(i,k) + phip(i,kp)*phi1(i,k))
          end do

        end if
      end do
    end do

    do i=1,np+2
      fac = fh(i)/rb(i)**2
      etheta(i,n) = etheta(i,n)*fac
      cnm(i,n)= cnm(i,n)*fac
      dnm(i,n)= dnm(i,n)*fac
    end do

  end do

  do i=1,np
    do n=1,nleg
      s1(i,n) = dnm(i,n) 
      s2(i,n) = cnm(i,n) + (etheta(i+1,n)-etheta(i-1,n))/(rb(i+1)-rb(i-1))
    end do
  end do
    
end subroutine hallterm

subroutine pol_matrix_inversion(dtb_myr,drb,eta_spe,aa,bb,cc)

  use grid, only: rb, np, belam, benu
  use legpol, only: nleg, frel
  use constants, only: PI
  use initial_magnetic, only: bcg
  implicit none
  integer :: i,n,k,kp,gg,j
  real*8 dtb_myr,drb(0:np+2),eta_spe(0:np+2,0:nleg)
  real*8 bcoef, icf, coef1, sum, fac1
  real*8 aa(nleg,nleg,np), bb(nleg,nleg,np), cc(nleg,nleg,np)

  do i=1,np
    fac1=(1d0-1d0/belam(i)**2)/(rb(i)*belam(i))
    do n=1,nleg    ! loop over n
      do k=1,nleg    ! loop over k
        sum = 0.d0
        do j=max(0,(n-nleg+k+1)/2),min(k,n)
          kp=k+n-2*j
          if (kp > 0) then
            gg = (kp+k+n)/2
            bcoef = bcg(gg-kp)*bcg(gg-k)*bcg(gg-n)/bcg(gg)
            ! Only axisymmetric terms (l=lp=m=0)
            icf=dsqrt((2*k+1)*(2*kp+1)*(2*n+1)/(4*PI))/(2*gg+1)*bcoef
            coef1 = 0.5d0*(n*(n+1)-kp*(kp+1)+k*(k+1))/dble(n*(n+1))
            sum = sum + coef1*icf*eta_spe(i,kp)
          end if
        end do      
    
        bb(n,k,i) = dtb_myr*benu(i)*sum*(2.d0/(belam(i)*drb(i))**2+k*(k+1)/rb(i)**2)
        aa(n,k,i) = - dtb_myr*benu(i)*sum*(1d0/(belam(i)*drb(i))**2 - fac1/drb(i) )
        cc(n,k,i) = - dtb_myr*benu(i)*sum*(1d0/(belam(i)*drb(i))**2 + fac1/drb(i) )
      end do      ! k loop
      
      coef1 = 1.d0/dsqrt(4.d0*PI)
      bb(n,n,i) = bb(n,n,i) +  &
     &          dtb_myr*benu(i)*coef1*eta_spe(i,0)*(2.d0/(belam(i)*drb(i))**2+n*(n+1)/rb(i)**2)
      aa(n,n,i) = aa(n,n,i) - dtb_myr*benu(i)*coef1*eta_spe(i,0)*(1d0/(belam(i)*drb(i))**2 - fac1/drb(i) )
      cc(n,n,i) = cc(n,n,i) - dtb_myr*benu(i)*coef1*eta_spe(i,0)*(1d0/(belam(i)*drb(i))**2 + fac1/drb(i) )

      bb(n,n,i) = 1.d0 + bb(n,n,i) 

    end do      ! n loop
  end do      ! i loop  (radial coordinate)


! boundary conditions 
  aa(:,:,1) = 0.d0
  do n=1,nleg
    do k=1,nleg
      bb(n,k,np) = bb(n,k,np) + cc(n,k,np)*(1.d0-k*frel(n)*drb(np)/rb(np+1))
    end do
  end do
  cc(:,:,np) = 0.d0

end subroutine pol_matrix_inversion


subroutine tor_matrix_inversion(dtb_myr,drb,eta_spe,aa,bb,cc)

  use grid, only: rb, np, belam, benu
  use legpol, only: nleg
  use constants, only: PI
  use initial_magnetic, only: bcg

  implicit none
  integer :: i,n,k,kp,gg,j
  real*8 dtb_myr,drb(0:np+2),eta_spe(0:np+2,0:nleg),deta_spe(0:np+1,0:nleg)
  real*8 bcoef, icf, coef1, sum, sumb, sumc
  real*8 aa(nleg,nleg,np), bb(nleg,nleg,np), cc(nleg,nleg,np)

  do i=1,np
    do n=0,nleg
      deta_spe(i,n)=(eta_spe(i+1,n)/belam(i+1)-eta_spe(i-1,n)/belam(i-1))/(2d0*drb(i)*belam(i))
    end do
  end do
  do i=1,np
    do n=1,nleg
      do k=1,nleg
        sum = 0.d0
        sumb = 0.d0
        sumc = 0.d0
        do j=max(0,(n-nleg+k+1)/2),min(k,n)
          kp=k+n-2*j
          iF (kp > 0) THEn
            gg = (kp+k+n)/2
            bcoef = bcg(gg-kp)*bcg(gg-k)*bcg(gg-n)/bcg(gg)
            icf=dsqrt((2*k+1)*(2*kp+1)*(2*n+1)/(4*PI))/(2*gg+1)*bcoef
            coef1 = 0.5d0*(n*(n+1)-kp*(kp+1)+k*(k+1))/dble(n*(n+1))
            sum = sum + coef1*icf*eta_spe(i,kp)
            sumb = sumb + coef1*icf*deta_spe(i,kp)
            sumc = sumc + icf*eta_spe(i,kp)
          endif
        end do      ! kp loop
  
        bb(n,k,i) = dtb_myr*benu(i)*( sum*(2.d0/(belam(i)*drb(i))**2 + sumc*k*(k+1)/rb(i)**2) )
        aa(n,k,i) = -dtb_myr*benu(i)*(sum/(belam(i)*drb(i))**2 - sumb/(2d0*drb(i)))
        cc(n,k,i) = -dtb_myr*benu(i)*(sum/(belam(i)*drb(i))**2 + sumb/(2d0*drb(i)))
      end do      ! k loop
     
      coef1 = 1.d0/dsqrt(4.d0*PI)
      bb(n,n,i) = bb(n,n,i) +   &
 &            dtb_myr*benu(i)*coef1*eta_spe(i,0)*(2.d0/(belam(i)*drb(i))**2+n*(n+1)/rb(i)**2)
      aa(n,n,i) = aa(n,n,i) -   &
 &            dtb_myr*benu(i)*coef1*(eta_spe(i,0)/(belam(i)*drb(i))**2 - deta_spe(i,0)/(2d0*drb(i)))
      cc(n,n,i) = cc(n,n,i) -   &
 &            dtb_myr*benu(i)*coef1*(eta_spe(i,0)/(belam(i)*drb(i))**2 + deta_spe(i,0)/(2d0*drb(i)))
      bb(n,n,i) = 1.d0 + bb(n,n,i) 
    end do      ! n loop
  end do      ! i loop  (radial coordinate)

  ! boundary conditions 

  aa(:,:,1) = 0.d0
  bb(:,:,1) = bb(:,:,1) + aa(:,:,1)   !  b.c. dPsi/dr=0
  cc(:,:,np) = 0.d0

end subroutine tor_matrix_inversion


subroutine triblocb(m,n,a,b,c,r,u)

  implicit none
  integer :: n,m
  real*8 :: a(n,n,m),b(n,n,m),c(n,n,m),r(n,m),u(n,m)
  integer :: i,j,k
  real*8 :: alfa(n,n),beta(n),alfinv(n,n)
  real*8 :: alf(n,n,m),bet(n,m)
  real*8 :: m1(n,n),m2(n,n),m3(n,n),m4(n,n)
  real*8 :: v1(n),v2(n),v3(n),v4(n),v5(n),v6(n)
  real*8 :: sum,dxmat(n,m)
  
  do j=1,n
    do i=1,n
      alf(i,j,1)=b(i,j,1)
    end do
    bet(j,1)=r(j,1)
  end do
  do k=2,m
    call transm(m,n,k-1,alf,alfa)
    call matinv(alfa,alfinv,n)  ! **   m1=a(k)
    call transm(m,n,k,a,m1) !  **   m2=c(k-1)
    call transm(m,n,k-1,c,m2) !  **   m3=a(k)*invalf(k-1)
    call matmul(m1,alfinv,m3,n,n,n) ! **   m4=a(k)*invalf(k-1)*c(k-1)
    call matmul(m3,m2,m4,n,n,n)       !  ** v1=a(k)*invalf(k-1)*beta(k-1)
    call transv(m,n,k-1,bet,beta)
    call mv_mul(m3,beta,v1,n)
    do i=1,n
    do j=1,n
        alf(i,j,k)=b(i,j,k)-m4(i,j)
      end do  
      bet(i,k)=r(i,k)-v1(i)
    end do  
  end do

  call transm(m,n,m,alf,alfa)
  call transv(m,n,m,bet,beta)
  call matinv(alfa,alfinv,n)
  call mv_mul(alfinv,beta,v1,n)   
  do i=1,n
    u(i,m)=v1(i)
  end do
  do k=m-1,1,-1
    call transm(m,n,k,alf,alfa)
    call transv(m,n,k,bet,beta)
    call matinv(alfa,alfinv,n)  !   ** m1=c(k)
    call transm(m,n,k,c,m1) !   ** v2=c(k)*u(k+1)
    call mv_mul(m1,v1,v2,n)  !  ** v3=beta(k)-c(k)*u(k+1)
    do i=1,n
       v3(i)=bet(i,k)-v2(i)
    end do  !     ** v1=alfinv(k)*( v3 )
    call mv_mul(alfinv,v3,v1,n)   !  ** u(k)=v1
    do i=1,n
      u(i,k)=v1(i) 
    end do
  end do  
  sum=0.
  do k=1,m
  call transm(m,n,k,a,m1)
  call transm(m,n,k,b,m2)
  call transm(m,n,k,c,m3)
  call transv(m,n,k-1,u,v1)
  call transv(m,n,k,u,v2)
  call transv(m,n,k+1,u,v3)
  call mv_mul(m1,v1,v4,n)
  call mv_mul(m2,v2,v5,n)
  call mv_mul(m3,v3,v6,n)
  do i=1,n
    dxmat(i,k)=v4(i)+v5(i)+v6(i)-r(i,k)
    sum=sum+dxmat(i,k)**2
  end do  
  end do
  sum=dsqrt(sum)

end subroutine triblocb
        

subroutine transm(m,n,k,a,b)
  implicit none
  integer :: i,j,k,m,n
  real*8 :: a(n,n,m),b(n,n)
  do i=1,n
  do j=1,n
     b(i,j)=a(i,j,k)
  end do
  end do   
end


subroutine transv(m,n,k,a,b)
  implicit none
  integer :: i,k,m,n
  real*8 :: a(n,m),b(n)
  do i=1,n
     b(i)=a(i,k)
  end do
end  

subroutine matinv(a,ainv,n)

  use math, only: lubksb, ludcmp
  implicit none
  integer :: n
  real*8 :: a(n,n), ainv(n,n)
  integer :: nphys, i, j
  parameter(nphys=200)
  real*8 :: temp(nphys,nphys), y(nphys,nphys), d
  integer :: index(nphys)

  ! make a copy of the array and initialize the identity matrix
  do j=1,n,1
    do i=1,n,1
      y(i,j) = 0.0
      temp(i,j) = a(i,j)
    end do
    y(j,j) = 1.0
  end do
  
  call ludcmp_b(temp,n,nphys,index,d)
  
  do j=1,n,1
    call lubksb_b(temp,n,nphys,index,y(1,j))
  end do
  
  do j=1,n,1
    do i=1,n,1
      ainv(i,j) = y(i,j)
    end do
  end do

end subroutine matinv
        
        
subroutine matmul(a,b,c,l,m,n)

  implicit none
  integer :: l, m, n
  real*8 :: a(l,m), b(m,n), c(l,n)
  integer :: i, j, k
  real*8 :: sum

  do i=1,l,1
    do j=1,n,1
      sum = 0.0
      do k=1,m,1
        sum = sum+a(i,k)*b(k,j)
      end do
      c(i,j) = sum
    end do
  end do

end subroutine matmul


subroutine mv_mul(a,v,rv,n)
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: a(n,n), v(n)
  real*8, intent(out) :: rv(n)
  integer :: i, j
  real*8 :: sum

  do i=1,n,1
    sum = 0.0
    do j=1,n,1
      sum = sum+a(i,j)*v(j)
    end do
    rv(i) = sum
  end do

end subroutine mv_mul


subroutine ludcmp_b(a,n,np,indx,d)
  implicit none
  integer :: n,np
  integer, dimension(n) :: indx
  real*8 :: d
  real*8, dimension(np,np) :: a
  real*8, parameter :: TINY=1.0d-25
  integer, parameter :: NMAX=1000
  integer :: i,imax,j,k
  real*8 :: aamax,dum,sum
  real*8, dimension(NMAX) :: vv

  d=1.d0
  do i=1,n
    aamax=0.d0
    do j=1,n
      if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
    end do
    if (aamax.eq.0.d0) then
      print*,'singular matrix in ludcmp'
      stop
    endif
    vv(i)=1.d0/aamax
  end do
  do j=1,n
    do i=1,j-1
      sum=a(i,j)
      do k=1,i-1
        sum=sum-a(i,k)*a(k,j)
      end do
      a(i,j)=sum
    end do 
    aamax=0.d0
    do i=j,n
      sum=a(i,j)
      do k=1,j-1
        sum=sum-a(i,k)*a(k,j)
      end do
      a(i,j)=sum
      dum=vv(i)*abs(sum)
      if (dum.ge.aamax) then
        imax=i
        aamax=dum
      endif
    end do
    if (j.ne.imax)then
      do k=1,n
        dum=a(imax,k)
        a(imax,k)=a(j,k)
        a(j,k)=dum
      end do
      d=-d
      vv(imax)=vv(j)
    endif
    indx(j)=imax
    if(a(j,j).eq.0.d0)a(j,j)=TINY
    if(j.ne.n)then
      dum=1.d0/a(j,j)
      do i=j+1,n
        a(i,j)=a(i,j)*dum
      end do
    endif
  end do
  
end subroutine ludcmp_b


subroutine lubksb_b(a,n,np,indx,b)
  implicit none
  integer :: n,np
  integer, dimension(n) :: indx
  real*8, dimension(np,np) :: a
  real*8, dimension(n) :: b
  integer :: i,ii,j,ll
  real*8 :: sum

  ii=0
  do i=1,n
    ll=indx(i)
    sum=b(ll)
    b(ll)=b(i)
    if (ii.ne.0)then
      do j=ii,i-1
        sum=sum-a(i,j)*b(j)
      end do
    else if (sum.ne.0.d0) then
      ii=i
    endif
    b(i)=sum
  end do
  do i=n,1,-1
    sum=b(i)
    do j=i+1,n
      sum=sum-a(i,j)*b(j)
    end do
    b(i)=sum/a(i,i)
  end do
  
end subroutine lubksb_b


subroutine interp(n,x,y,x0,y0)
  implicit none
  integer j,n,jl,jm,ju
  real*8 x(n), y(n), x0, y0, alfa
  if (x0 >= x(n)) then
    y0 = y(n)
    return
  endif
  jl=0
  ju=n+1

  do while (ju-jl > 1)
    jm=(ju+jl)/2
    if((x(n).gt.x(1)).eqv.(x0.gt.x(jm)))then
      jl=jm
    else
      ju=jm
    endif
  end do

  j=jl
  if ( j == 0 .or. j == n ) then
     write(*,*)'out of the table', j, n, x0
     stop
  end if
  alfa = (x0-x(j))/(x(j+1)-x(j))
  y0 = y(j) + alfa*(y(j+1)-y(j))

end subroutine interp
  
