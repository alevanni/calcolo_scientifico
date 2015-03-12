function II=inversa(diago, ssd, n)
%calcola l'inversa di una matrice tridiagonale simmetrica
%con lo stesso valore lungo tutta la diagonale
II=zeros(n,n);
ei=zeros(n,1); %vettore della base canonica
z=zeros(n,1);
x=zeros(n,1);
%si calcola la fattorizzazione LU (non riportata)
[l,u,v]=LU(diago, ssd, n);
for i=1:n    !si risolvono i sistemi con la base canonica come termine noto
%si inizializza il vettore e_i 
ei(i)=1;
z=lowert(l, ei, n );
x=uppert(u,v,z,n);
II(:,i)=x;
%si torna al vettore di soli zeri
ei(i)=0;
endfor

endfunction
