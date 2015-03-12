subroutine gradientep2( n,  iter, sol, b)
  ! esegue iter passi dell'iterazione del gradiente coniugato
  ! per risolvere il sistema Ax=b con A matrice del laplaciano discreto
  ! precondizionato con la matrice P2 (diagonale a blocchi tridiagonali)
  implicit none
  integer, parameter :: dp=KIND(0.d0)
  integer(16) :: n, iter !dimensione dei blocchi e del sistema
 
  real(dp), dimension(1:iter) :: r, b, x, p, q, sol, diff

  real(dp), dimension(1:iter) :: z!soluzione di P2*z=r
  real(dp), dimension(1:iter) :: aux !soluzione intermedia di P1*z=r
  integer(16) ::  i, j
  real(dp), dimension(1:n) ::  u    !
  real(dp) :: dis, diag !rispettivamente numero sulla diagonale inferiore 
  !(e superiore) di un blocco di P2 e numero sulla diagonale
  real(dp), dimension(1:n-1) :: l,  v !rispettivamente sottodiagonale
  ! di un blocco di L , diagonale e sovradiagonale di U
  real(dp) :: rho, alpha, beta, rhop,  nor1, nor
   
  !definisco un blocco di P2

  diag=2
  dis=-1
  
  ! Calcolo il residuo r = b - A x
  ! che per x=0 coincide con b
 
  x=0
  r=b
  call LU(diag, dis, n, l, u, v) !calcolo una sola volta la fattorizzazione
  ! LU di un blocco
  do i=1,iter

     !risolvo il sistema P2*z=r a blocchi di dimensione n
          
     do j=0,n-1  
     !risolvo L*aux=r
     call lower(l, r(j*n+1:(j+1)*n), n, aux(j*n+1:(j+1)*n)) 
     !risolvo U*z=aux
     call upper(u, v, aux(j*n+1:(j+1)*n), n, z(j*n+1:(j+1)*n) ) 
     end do

     rho = SUM(r*z) !norma quadra del resto

     if(i==1) then
        p = z
     else
        beta = rho / rhop 
        p = z + beta * p
     end if

     
     ! Calcolo q = Ap
     call Ax(n, p, q) 
        
     alpha = rho / SUM(p*q) !prodottoscalare

     ! aggiorno le variabili per l'iterazione successiva
     x = x + alpha * p

     ! Calcolo il prossimo residuo; r = b - Ax ; ora x = x + alpha * p
     ! e quindi r diventa b - A(x + alpha*p) = b - Ax - alpha*q =
     ! r - alpha * q

     r = r - alpha * q
     rhop = rho
     diff=x-sol !differenza tra la soluzione reale e quella calcolata
     !devo calcolare || x-r||
     nor1 = sqrt(sum(r*r))
     nor=sqrt(sum(diff*diff))
     open(unit=8, file="residuoP2.txt")
     write(8, *) nor1
     open(unit=7, file="differenzaP2.txt")
     write(7, *) nor
     !controllo il residuo
     if( nor < 1.D-7 ) then
        exit
     end if
  end do
end subroutine gradientep2
