subroutine lower(l, y, n, z)
!risolve un sistema bidiagonale superiore. In entrata prende i valori della 
!sottodiagonale (la diagonale e' composta da soli 1), il termine noto e la
! dimensione del sistema;
implicit none 
integer, parameter :: dp=kind(0.d0)
integer(16):: n, i !dimensione, indice
real(dp), dimension (1:n) :: z, y !rispettivamente soluzione (da trovare) e 
!termine noto
real(dp), dimension (1:n-1) :: l !sottodiagonale
!si inizia dalla cima
z(1)=y(1)

do i=2,n
z(i)=y(i)-l(i-1)*z(i-1)
end do


end subroutine lower
