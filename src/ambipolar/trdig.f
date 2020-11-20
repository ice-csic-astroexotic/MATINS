c    This program generates the block-tridiagonal system
c
c    | b1 c1 0  0  0.....| |x1|    |r1|
c    | a2 b2 c2 0  0.....| |x2|    |r2|
c    | 0  a3 b3 c3 0.....| |. |  _ |. |
c    | ..................| |. |  _ |. |
c    | ..................| |. |    |. |
c    | .............an bn| |xn|    |. |
c
c

      subroutine trdig(nx,nz,pci,cc,s)
            
      use math, only: solvtb

      implicit none
      integer m,n,nx,nz
      integer j1,j2,j3,ik1,ik2
      real*8 a(nx,nz,3),b(nx,nz,3),c(nx,nz,3)
      real*8 x(nx,nz),r(nx,nz)
      real*8 cc(nx,nz,5)
      real*8 pci(nx,nz),s(nx,nz)             

c     initialize a,b,c,x,b to zero

         x=0.d0
         r=0.d0
         a=0.d0
         b=0.d0
         c=0.d0

c  hacemos el cambio 
c                     
c       m,n-1   ----->  1
c       m-1,n   ----->  2  
c       m,n     ----->  3
c       m+1,n   ----->  4 
c       m,n+1   ----->  5 

      do m=1,nx
      do n=1,nz          
         a(m,n,2)=cc(m,n,2)
         b(m,n,1)=cc(m,n,1)
         b(m,n,2)=cc(m,n,3)
         b(m,n,3)=cc(m,n,5)
         c(m,n,2)=cc(m,n,4)
         r(m,n)=s(m,n)
      enddo
      enddo
c      write(*,*)'trdlg'      
c      write(*,*)'a',((a(m,n,1),m=1,nx),n=1,nz)
c      write(*,*)'c',((c(m,n,1),m=1,nx),n=1,nz)
c      write(*,*)'b1',((b(m,n,1),m=1,nx),n=1,nz)
c      write(*,*)'b2',((b(m,n,2),m=1,nx),n=1,nz)
c      write(*,*)'b3',((b(m,n,3),m=1,nx),n=1,nz)

      call solvtb(nx,nz,a,b,c,r,x)

123   format(5(1x,d10.4))     
      do m=1,nx
      do n=1,nz
         pci(m,n)=x(m,n)
       enddo
       enddo
       

      return
      end
        
       
