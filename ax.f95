subroutine Ax(n, x, y)
! Calcola il prodotto tra A (matrice del laplaciano discreto) e un vettore, 
! che viene memorizzato in b. Prende in input la dimensione dei blocchi di A
  implicit none
  integer, parameter:: dp=KIND(0.d0)
  integer(16) :: n, i,j
  
  real(dp), dimension(1:n*n)   :: x, y
  real(dp), dimension(1:n , 1:n) :: xm, ym !il prodotto verra' 
  !eseguito a blocchi
  
 
  xm=reshape(x, (/n, n/) ) !adesso x e' una matrice nxn
  ym=reshape(y, (/n, n/) ) !adesso y e' una matrice nxn
   
  
  do i=1,n
      
     do j=1,n
         
      if (i==1) then !prima riga
         
          if (j==1) then
            ym(i,j)=4*xm(i,j) -xm(i+1,j) -xm(i,j+1)
          else if (j==n) then 
            ym(i,j)=4*xm(i,j)  -xm(i+1,j) - xm(i,j-1)
          else 
          ym(i,j)=4*xm(i,j) -xm(i+1,j) - xm(i,j-1)-xm(i,j+1)
          end if     
     


     else if (i==n) then  !ultima riga
        
        if (j==1) then
          ym(i,j)=4*xm(i,j) - xm(i-1,j)  -xm(i,j+1)
        else if (j==n) then 
          ym(i,j)=4*xm(i,j) - xm(i-1,j)  - xm(i,j-1)
        else 
         ym(i,j)=4*xm(i,j) - xm(i-1,j)  - xm(i,j-1)-xm(i,j+1)
       end if
       
     else if (j==n .and. i>1) then !ultima colonna
         ym(i,j)=4*xm(i,j) - xm(i-1,j) -xm(i+1,j) - xm(i,j-1)
     else if (j==1 .and. i>1) then !prima colonna
         ym(i,j)=4*xm(i,j) - xm(i-1,j) -xm(i+1,j) -xm(i,j+1)


     else
      ym(i,j)=4*xm(i,j) - xm(i-1,j) -xm(i+1,j) - xm(i,j-1)-xm(i,j+1)
     end if  
     end do
   end do
 x=reshape(xm, (/n*n/)) !si torna ai vettori
 y=reshape(ym, (/n*n/))
end subroutine ax
