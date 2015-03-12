subroutine upper(u, v, z, n, x)
!risolve un sistema bidiagonale superiore. In entrata prende i vettori
! della diagonale (u),
!della sopradiagonale (v), il termine noto (z) e la dimensione (n)
!memorizza la soluzione nel vettore x
implicit none 
integer, parameter :: dp=kind(0.d0)
integer(16):: n, i !dimensione, indice
real(dp), dimension (1:n) :: x, z, u !soluzione, termine noto e diagonale
real(dp), dimension (1:n-1) :: v !sopradiagonale

!si inizia dal fondo
x(n)=z(n)/u(n)

do i=n-1,1,-1
x(i)=( z(i)-v(i)*x(i+1) )/u(i)
end do
end subroutine upper
