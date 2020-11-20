IF (bcout .eq. 0) then
  CALL VACUUMBC()

  do k=2,kmax
    i=2*k-2
    thsup(k-1)=theta(i)
    bsup(k-1)=br(i,np)
        enddo
  
     ELSEIF (bcout .eq. 1) then
 CALL INVERT_BC(kmax-1,bsup,btheta)
 do k=1,kmax
  i=2*k-1
  bth(i,np+1)=btheta(k)
 enddo



 ; Initial:
SUBROUTINE GET_GREEN(m)
  implicit none
  integer i,j,k,m,mi,resi
  parameter (resi=10)
  real*8 pi
  real*8 dthm,thm(m)
  real*8 fw(1000,1000),fint
  real*8 thi(resi*m),dthi
  common /fw_green/ fw,dthm

  mi=resi*m
  if (mi .gt. 1000) then
write(6,*) 'GET_GREEN: resi too high!'
stop
  endif
  pi=dacos(-1.d0)
  dthm=pi/dble(m)
  dthi=pi/dble(mi)
  do j=1,m
thm(j)=dble(j-0.5d0)*dthm
  enddo
  do k=1,mi
thi(k)=dble(k-0.5d0)*dthi
  enddo

  fw=0.d0
  do i=1,m
  do j=1,m
  do k=1,mi
if (thi(k) .ge. thm(j)-0.5d0*dthm .and.
 &	    thi(k) .le. thm(j)+0.5d0*dthm) then
fw(i,j)=fw(i,j)+fint(thm(i),thi(k),pi)*dthi
endif
  enddo
  enddo
  enddo

  return
  END


  SUBROUTINE INVERT_BC(m,brad,btheta)
    implicit none
    integer i,j,m
    integer*4 indx(m)
    real*8 pi,d
    real*8 brad(m)
    real*8 dthm,btheta(1:m+1)
    real*8 matpsi(m,m),b(m)
    real*8 fw(1000,1000)
    common /fw_green/ fw,dthm

    pi=dacos(-1.d0)
    do i=1,m
    do j=1,m
matpsi(i,j)=0.5d0*fw(i,j)
if (j .eq. i)
   &	matpsi(i,j)=matpsi(i,j)+2.d0*pi
    enddo
    enddo

    b=0.d0
    do i=1,m
    do j=1,m
      b(i)=b(i)-brad(j)*fw(i,j)
    enddo
    enddo

    call ludcmp(matpsi,m,m,indx,d)
    call lubksb(matpsi,m,m,indx,b)

    do i=2,m
btheta(i)=(b(i)-b(i-1))/dthm
    enddo
    btheta(1)=0.d0
    btheta(m+1)=0.d0

    return
    END
