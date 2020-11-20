!-----------------------------------------------------------------------
! Set up the operator and source, define the boundary conditions, and
! solve the linear system.
!
! Called by:
! FORCE_FREE
!
! Calls:
! MULTIPOLE (Multipole expansion.)
! TRDIG (Solution of linear system.)
!-----------------------------------------------------------------------
	subroutine OPERATOR(Pold,Pnew,rd,td,cth,sth,dr,dt,geo,nx,nz,
     &	factor,zstar,corr,dpl,lmax,Pc,at,ma,fit_type)

	implicit none

! Dimension indices.
	integer m,n,nx,nz

! Radial and angular arrays, and related coordinates.
	real*8 rd(nx),td(nz),cth(nz),sth(nz)
	real*8 dr,dt

! Operator, source and poloidal function.
	real*8 a(nx,nz,5),Q(nx,nz),Q4(nz),Q5(nz)
	real*8 Pold(nx,nz),Pnew(nx,nz)
	real*8 Pave(nz)

! Number of multipoles and related variables.
	integer l,lmax
	real*8 al(lmax),dpl(lmax,nz)

! Toroidal function.
	integer fit_type,ma
	real*8 at(ma),dft,ft,x
	real*8 Pc

! Average correction wrt the guess.
	real*8 corr

! Relativistic case.
	real*8 z,zstar
	real*8 fn

! Source (in the Grad-Shafranov equation).
	real*8 factor

! Geometric factor for an uneven grid.
	real*8 drm,geo

!-----------------------------------------------------------------------
! General matrix elements (Grad-Shafranov equation).
	a=0d0
	do n=2,nz-1
	   do m=2,nx-1
! Gometric grid step.
	      drm=dr*geo**(m-2)
! Relativistic correction.
	      z=zstar/rd(m)
! Coefficient of P(m,n-1).
	      a(m,n,1)=1d0/(rd(m)*dt)**2
     &		 +cth(n)/(2d0*rd(m)**2*sth(n)*dt)
! Coefficient of P(m-1,n).
	      a(m,n,2)=2d0/((geo+1d0)*drm**2)
cc	      a(m,n,2)=(1d0-z)/dr**2-z/(2d0*rd(m)*dr)
! Coefficient of P(m,n).
! Choice 1: Toroidal field as an operator.
! Choice 1.1: Linear function everywhere, T=sP.
cc	      a(m,n,3)=-2d0/dr**2-2d0/(rd(m)*dt)**2+s**2
! Choice 1.2: Step function, T=s(P-Pc) for P>Pc, and zero otherwise.
cc	      if(Pold(m,n).gt.Pc)then
cc		 a(m,n,3)=-2d0/dr**2-2d0/(rd(m)*dt)**2+s**2
cc	      else
cc		 a(m,n,3)=-2d0/dr**2-2d0/(rd(m)*dt)**2
cc	      endif
! Choice 2: Toroidal field as a source. (Assign the appropriate source.)
	      a(m,n,3)=-2d0/(geo*drm**2)-2d0/(rd(m)*dt)**2
cc	      a(m,n,3)=-2d0*(1d0-z)/dr**2-2d0/(rd(m)*dt)**2
! Coefficient of P(m+1,n).
	      a(m,n,4)=2d0/(geo*(geo+1d0)*drm**2)
cc	      a(m,n,4)=(1d0-z)/dr**2+z/(2d0*rd(m)*dr)
! Coefficient of P(m,n+1).
	      a(m,n,5)=1d0/(rd(m)*dt)**2
     &		 -cth(n)/(2d0*rd(m)**2*sth(n)*dt)
	   enddo
	enddo

!-----------------------------------------------------------------------
! Boundary conditions.
! Axial boundary conditions (n=1 and n=nz).
	do m=1,nx
! North pole.
	   a(m,1,1)=0d0
	   a(m,1,2)=0d0
	   a(m,1,3)=1d0
	   a(m,1,4)=0d0
	   a(m,1,5)=0d0
! South pole.
	   a(m,nz,1)=0d0
	   a(m,nz,2)=0d0
	   a(m,nz,3)=1d0
	   a(m,nz,4)=0d0
	   a(m,nz,5)=0d0
	enddo
! Radial boundary conditions.
	do n=2,nz-1
! Inner boundary (m=1).
	   a(1,n,1)=0d0
	   a(1,n,2)=0d0
	   a(1,n,3)=1d0
	   a(1,n,4)=0d0
	   a(1,n,5)=0d0
! Outer boundary (m=nx).
! Dirichlet boundary conditions (for the function).
! Type 1: Specify the value of the function at the boundary.
	   a(nx,n,1)=0d0
	   a(nx,n,2)=0d0
	   a(nx,n,3)=1d0
	   a(nx,n,4)=0d0
	   a(nx,n,5)=0d0
	enddo

! Alternative boundary conditions for the outer boundary.
! Neumann boundary conditions (for the derivatives).
! Type 2: Neumann boundary conditions of the first kind.
! Rewrite the operator at the boundary, and set the source to zero.
! Outer boundary (m=nx).
cc	do n=2,nz-1
cc	   a(nx,n,1)=1d0/(rd(nx)*dt)**2
cc     &	      +cth(n)/(2d0*rd(nx)**2*sth(n)*dt)
cc	   a(nx,n,2)=2d0/dr**2
cc	   a(nx,n,3)=-2d0/dr**2-2d0/(rd(nx)*dt)**2-2d0*ell/(rd(nx)*dr)
cc	   a(nx,n,4)=0d0
cc	   a(nx,n,5)=1d0/(rd(nx)*dt)**2
cc     &	      -cth(n)/(2d0*rd(nx)**2*sth(n)*dt)
cc	enddo
! Type 3: Neumann boundary conditions of the second kind.
! Rewrite the operator at the boundary, and update the source.
! Outer boundary (m=nx).
cc	do n=2,nz-1
cc	   a(nx,n,1)=1d0/(rd(nx)*dt)**2
cc     &	      +cth(n)/(2d0*rd(nx)**2*sth(n)*dt)
cc	   a(nx,n,2)=2d0/dr**2
cc	   a(nx,n,3)=-2d0/dr**2-2d0/(rd(nx)*dt)**2
cc	   a(nx,n,4)=0d0
cc	   a(nx,n,5)=1d0/(rd(nx)*dt)**2
cc     &	      -cth(n)/(2d0*rd(nx)**2*sth(n)*dt)
cc	enddo
! Type 4: Neumann boundary conditions of the third kind.
! Impose numerical derivative at the boundary using the previous point.
! The operator remains diagonal.
cc
! Type 5: Neumann boundary conditions of the fourth kind.
! Impose numerical derivative at the boundary using the previous point.
! Similar to the type 4, but with some terms moved to the operator.
	drm=dr*geo**(nx-2)
	do n=2,nz-1
	   a(nx,n,1)=0d0
	   a(nx,n,2)=-1d0/drm
	   a(nx,n,3)=1d0/drm
	   a(nx,n,4)=0d0
	   a(nx,n,5)=0d0
	enddo

!-----------------------------------------------------------------------
! Multipole expansion at m=nx-1.
! Q4 needs to be calculated at nx-1, while Q5 can be calculated at
! either nx-1 or nx (to the same accuracy), or as some average (better).
! If P decreases as 1/r (for a dipole), the first order derivative at
! the outer boundary is exactly equal to the geometric average of P_j-1
! and P_j at the geometric average of r_j-1 and r_j.
	Pave = sqrt(Pold(nx-1,:)*Pold(nx,:))
	call MULTIPOLE(Pave,cth,dpl,lmax,nz,al)
cc	call MULTIPOLE(Pold(nx-1,:),cth,dpl,lmax,nz,al)

cc	z=zstar/rd(nx-1)
cc	open(1,file="out/frel.txt")
cc	do l=1,lmax
cc	   call GET_FREL(l,z,fn)
cc	   write(1,'(i5,e17.8)')l,fn
cc	enddo
cc	close(1)
cc	stop

! Boundary values for type 4 and type 5.
	Q4=0d0
	Q5=0d0
! Relativistic factors (z and fn).
	z=zstar/rd(nx-1)
	do l=1,lmax
	   call GET_FREL(l,z,fn)
	   do n=1,nz
cc	      Q4(n)=Q4(n)+(1d0-cth(n)**2)
cc     &		 *al(l)*(rd(nx-1)/rd(nx))**l*dpl(l,n)
cc	      Q5(n)=Q5(n)-(1d0-cth(n)**2)
cc     &		 *al(l)*(dble(l)/rd(nx-1))*dpl(l,n)
cc	      Q5(n)=Q5(n)-(1d0-cth(n)**2)
cc     &		 *fn*al(l)*(dble(l)/rd(nx-1))*dpl(l,n)
! Geometric average (for a dipole).
	      Q5(n)=Q5(n)-(1d0-cth(n)**2)
     &		 *fn*al(l)*(dble(l)/sqrt(rd(nx-1)*rd(nx)))*dpl(l,n)
	   enddo
	enddo

!-----------------------------------------------------------------------
! Source.
! Depends how the operator is defined and on the boundary conditions!
	Q=0d0

! Toroidal field.
	do m=2,nx-1
	   do n=2,nz-1
! Choice 1: Toroidal field as an operator.
! Choice 1.1: Linear function everywhere, T=sP.
cc	      Q(m,n)=0d0
! Choice 1.2: Step function, T=s(P-Pc) for P>Pc, and zero otherwise.
cc	      if(Pold(m,n).gt.Pc)then
cc		 Q(m,n)=s**2*Pc
cc	      endif
! Choice 2: Toroidal field as a source.
! Choice 2.1: Linear function everywhere, T=sP.
cc	      Q(m,n)=-s**2*Pold(m,n)
! Choice 2.2: Step function, T=s(P-Pc)^sigma for P>Pc, and zero
! otherwise.
cc	      if(Pold(m,n).gt.Pc)then
cc		 Q(m,n)=-s**2*sigma*(Pold(m,n)-Pc)**(2d0*sigma-1d0)
cc		 Q(m,n)=-s**2*(Pold(m,n)-Pc)
cc	      endif
! Toroidal function as an external function.
	      x=Pold(m,n)
	      if(x.gt.Pc)then
	         call FTOR(x,at,ft,dft,ma,fit_type)
	         Q(m,n)=-ft*dft
	      endif
! Density term (rho = 1 - r^2) as a source (for the interior).
cc	      if(rd(m).le.1d0)then
cc		 Q(m,n)=Q(m,n)+factor*(1d0-rd(m)**2)*rd(m)**2*sth(n)**2
cc	      endif
	   enddo
! Relativistic correction.
	   z=zstar/rd(m)
	   Q(m,:)=Q(m,:)/(1d0-z)
	enddo

!-----------------------------------------------------------------------
! Axial boundaries (at n=1 and n=nz).
! The poloidal function along the axis is zero.
	Q(:,1)=0d0
	Q(:,nz)=0d0

! Radial boundaries (at m=1 and m=nx).
! Inner boundary.
	Q(1,:)=Pold(1,:)
cc	Q(1,:)=0d0

! Outer boundary.
! Type 1: Dirichlet boundary conditions (for the function).
cc	Q(nx,:)=Pold(nx,:)
cc	Q(nx,:)=0d0
! Type 2: Neumann boundary conditions of the first kind.
cc	Q(nx,:)=0d0
! Type 3: Neumann boundary conditions of the second kind.
cc	Q(nx,:)=2d0*dble(ell)*Pold(nx,:)/(rd(nx)*dr)
! Type 4: Neumann boundary conditions of the third kind.
! Impose numerical derivatives at the boundary using the inner points.
! The derivative is added as a source.
! Single multipole.
cc	Q(nx,:)=Pold(nx-1,:)*(rd(nx-1)/rd(nx))**ell
! Multipole expansion.
cc	Q(nx,:)=Q4(:)
! Type 5: Neumann boundary conditions of the fourth kind.
	Q(nx,:)=Q5(:)

!-----------------------------------------------------------------------
! Solution of the linear system.
	Pnew=0d0
	call TRDIG(nx,nz,Pnew,a,Q)

! Calculate the average correction wrt the guess.
	corr=0d0
	do m=1,nx
	   do n=1,nz
	      corr=corr+(Pnew(m,n)-Pold(m,n))**2
	   enddo
	enddo
	corr=corr/dble(nx*nz)

	return
	end
