subroutine Ax2(n, x, y)
    ! Calcola cioe' il prodotto tra A (matrice del laplaciano discreto)
    ! e un vettore,
! che viene memorizzato in b
  implicit none
  integer, parameter :: dp=KIND(0.d0)
  integer(16) :: n, i, j
  
  real(dp),dimension(1:n*n)   :: x, y
  real(dp),dimension(1:n , 1:n) :: xm, ym !il prodotto verra' 
  !eseguito a blocchi
  
  
 
  xm=reshape(x, (/n, n/) )!adesso x e' una matrice nxn
  ym=reshape(y, (/n, n/) ) !adesso y e' una matrice nxn

  !angoli
  ym(1,1)=4*xm(1,1) -xm(2,1) -xm(1,2)
  ym(1,n)=4*xm(1,n)  -xm(2,n) - xm(1,n-1)
  ym(n,1)=4*xm(n,1) - xm(n-1,2)  -xm(n,2)
  ym(n,n)=4*xm(n,n) - xm(n-1,n)  - xm(n,n-1)  
  
!prima riga e ultima riga insieme
do j=2,n-1
ym(1,j)=4*xm(1,j) -xm(2,j) - xm(1,j-1)-xm(1,j+1)
ym(n,j)=4*xm(n,j) - xm(n-1,j)  - xm(n,j-1)-xm(n,j+1)
end do

!prima colonna e ultima colonna insieme
do i=2,n-1
ym(i,1)=4*xm(i,1) - xm(i-1,1) -xm(i+1,1) -xm(i,2)
ym(i,n)=4*xm(i,n) - xm(i-1,n) -xm(i+1,n) - xm(i,n-1)
end do
!sarebbe stato ancora piu' efficiente se avessi accorpato i cicli precedenti!

!tutto il resto

do i=2,n-1
do j=2,n-1
ym(i,j)=4*xm(i,j) - xm(i-1,j) -xm(i+1,j) - xm(i,j-1)-xm(i,j+1)
end do
end do


  
 x=reshape(xm, (/n*n/))
 y=reshape(ym, (/n*n/))
end subroutine Ax2
