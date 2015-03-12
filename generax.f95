subroutine generax(n, x)  !genera a caso un vettore di dimensione n
implicit none
integer, parameter :: dp=KIND(0.d0)

integer(16) :: n !dimensione
real(dp), dimension (1:n) :: x !soluzione


call random_number(x)
end subroutine
