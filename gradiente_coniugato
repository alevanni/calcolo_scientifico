subroutine gradiente_coniugato( n,  iter, sol, b)
  ! esegue iter passi dell'iterazione di gradiente coniugato
  ! per risolvere il sistema Ax=b con A matrice del laplaciano discreto
  ! In entrata passo la dimensione dei blocchi (serve per il prodotto 
  ! matrice-vettore), quella della matrice, la soluzione reale e il termine
  ! noto
  implicit none
  integer, parameter :: dp=KIND(0.d0)
  integer(16) :: n, i, iter 
  real(dp), dimension(1:iter) :: r, b, x, p, q, sol, diff  
  real(dp) :: rho, alpha, beta, rhop,  nor1, nor
 
  ! Calcolo il residuo r = b - A x
  ! che per x=0 coincide con b
  x=0
  r=b
  do i=1,iter
     rho = SUM(r*r) !norma quadra del resto
     if(i==1) then
        p = r
     else
        beta = rho / rhop 
        p = r + beta * p
     end if

     
     ! Calcolo q = Ap
     call Ax(n, p, q) 

     alpha = rho / SUM(p*q) !prodottoscalare

     ! aggiorno le variabili per passo l'iterazione successiva
     x = x + alpha * p

     ! Calcolo il prossimo residuo; r = b - Ax ; ora x = x + alpha * p
     ! e quindi r diventa b - A(x + alpha*p) = b - Ax - alpha*q =
     ! r - alpha * q
     r = r - alpha * q
     rhop = rho
     diff=x-sol !differenza tra la soluzione real e quella calcolata
     !devo calcolare || x-r||
     
     nor1 = sqrt(sum(r*r))
     nor=sqrt(sum(diff*diff))
     open(unit=1, file="residuo.txt")
     write(1, *) nor1
     open(unit=3, file="differenza.txt")
     write(3, *) nor
     !controllo il residuo
     if( nor < 1.D-7 ) THEN
     ! write(*, *) i, nor1  scommentare per stampare il 
     ! minimo numero di passi
           exit
     end if
  end do
end subroutine gradiente_coniugato
