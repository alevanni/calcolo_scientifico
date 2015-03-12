subroutine gradientep1(n, iter, sol, b)
  ! esegue iter passi dell'iterazione del gradiente coniugato
  ! per risolvere il sistema Ax=b con A matrice del laplaciano discreto
  ! precondizionato con la matrice P1 (n^2xn^2, con 4 sulla diagonale e -1
  ! sulla sopra/sottodiagonale)
  
  implicit none
  integer, parameter :: dp=kind(0.d0)
  integer(16) :: n, iter, i
  real(dp) :: dis, diag ! rispettivamente numero sulla diagonale
  ! inferiore (e superiore) 
  !e numero sulla diagonale
  real(dp), dimension(1:iter) :: r, b, x, p, q, sol, diff
  real(dp), dimension(1:iter) :: z !soluzione di P1*z=r
  real(dp), dimension(1:iter) :: aux !soluzione intermedia di P1*z=r
  real(dp), dimension(1:iter) ::  u !diagonale di U 
  real(dp), dimension(1:iter-1) :: l,  v !rispettivamente sottodiagonale 
  ! di L e sovradiagonale di U
  real(dp) :: rho, alpha, beta, rhop,  nor1, nor
   
  !definisco P1 

  diag=4
  dis=-1
  
  ! Calcolo il residuo r = b - A x
  ! che per x=0 coincide con b
  
  x=0
  r=b
  !calcolo la fattorizzazione LU di P1
  call LU(diag, dis, iter, l, u, v) 
  do i=1,iter

     !risolvo il sistema P1*z=r 
     
     call lower(l, r, iter, aux)     !risolvo L*aux=r
     call upper(u, v, aux, iter, z ) !risolvo U*z=aux. 
     ! Adesso ho trovato z.
     rho = SUM(r*z) 
     if(i==1) then
        p = z
     else
        beta = rho / rhop 
        p = z + beta * p
     end if

     
     ! Calcolo q = Ap
     call Ax(n, p, q) 

     alpha = rho / SUM(p*q) !prodottoscalare

     !aggiorno le variabili per passo l'iterazione successiva
     x = x + alpha * p

     ! Calcolo il prossimo residuo; r = b - Ax ; ora x = x + alpha * p
     ! e quindi r diventa b - A(x + alpha*p) = b - Ax - alpha*q =
     ! r - alpha * q
     r = r - alpha * q
     rhop = rho
     diff=x-sol !differenza tra la soluzione reale e quella calcolata
     !devo calcolare || x-r||
     ! Controllo il residuo
     nor1 = sqrt(sum(r*r))
     !Controllo la differenza
     nor=sqrt(sum(diff*diff))
      open(unit=4, file="residuoP1.txt")
     write(4, *) nor1
     open(unit=5, file="differenzaP1.txt")
     write(5, *) nor
     if( nor < 1.D-7 ) then        
        exit
     end if
  end do
  
end subroutine gradientep1
