subroutine LU(d, dis, n, l, u, v) 
!calcola la decomposizione LU di una matrice tridiagonale simmetrica,
! con lo stesso valore lungo tutta la diagonale e sovra/sottodiagonale; 
!prende in input due numeri (numero sulla diagonale e numero sulla 
!diagonale inferiore) e la dimensione della matrice
!restituisce tre vettori: l (sottodiagonale di L), u e v 
!(rispettivamente diagonale e sovradiagonale di U)

implicit none

integer, parameter :: dp=kind(0.d0)
integer(16):: n, i !dimensione, indice
real(dp) :: dis, d !elemento sulla diagonale inferiore (e superiore) 
!e sulla diagonale
real(dp), dimension (1:n) ::  u !diagonale di u
real(dp), dimension (1:n-1) ::  l, v !sottodiagonale di L e 
!sovradiagonale di U

!primo passo 

u(1)=d
v(1)=dis
l(1)=dis/u(1)

do i=2,n-1
u(i)=d-l(i-1)*v(i-1)
v(i)=dis
l(i)=dis/u(i)
end do
!ultimo passo
u(n)=d-l(n-1)*v(n-1)
!write(*,*) "ok"
end subroutine LU
