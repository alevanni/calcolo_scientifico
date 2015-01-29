program main1
implicit none

INTEGER, PARAMETER :: dp=KIND(0.d0)
INTEGER(16) ::  n, iter
real(dp), allocatable :: sol(:), b(:), x(:) 

!chiede la dimensione
write(*,*) "Dammi la dimensione dei blocchi"
read(*,*) n

iter=n*n

allocate(x(iter), b(iter), sol(iter)) 
call generax(iter,sol) !genero a caso la soluzione
call Ax(n,sol,b) !trovo il termine noto b
!call gradiente_coniugato(n, iter, sol, b) !cerco la soluzione x del sistema Ax=b
call gradientep(n, iter, sol, b)

end program main1
