function [U Bi V] = Bidiag_Francis_Step_Update_U_V(U, V, Bi)
%copy Bidaig_Francis_Step here
%Francis step to chase the bulge to diagonalize matrix called from implicit_bidaig_QR_SVD
%10.3.5.1
[m, n] = size(Bi);
if(m~=n)  
   return; %this code works for initial call from implicit_bidiag_qr_svd. then odd matrix
end
assert( m~=n, 'only square diagonal matrix taken as input' );
%U = eye(m); %accunulate multiplied val of Go*G1...Take from input/multiply
%V =eye(m);%accunulate multiplied values of Go'*G1'...
%G1 w T's elements 
t00 = Bi(1,1) * Bi(1,1);
t10 = Bi(1,2) * Bi(1,1);
tm1m1 = ( Bi(m-1,m) * Bi(m-1,m) ) * ( Bi(m,m) * Bi(m,m));
%tm1m1 = ( Bi(m-1,n) * Bi(m-1,n) ) * ( Bi(m,n) * Bi(m,n));
G = eye(2);
IG = eye (m);
for j = 1:m-1
    if j == 1
        G = Givens_rotation( [ t00 - tm1m1;t10 ]);  
        U([j,j+1],:) = G'*U([j, j+1],:);    
        V(:,[j j+1]) = V(:,[j j+1])*G;
     elseif (j == m-1)
        G= Givens_rotation([Bi(j,j);Bi(j+1,j)]);
        U([j-1,j],:) = G'*U([j-1, j],:);     
        V(:,[j-1 j]) = V(:,[j-1 j])*G;
    else 
        G= Givens_rotation([Bi(j,j);Bi(j+1,j)]);
        U([j-1,j+1],:) = G'*U([j-1, j+1],:);    
        V(:,[j-1 j+1]) = V(:,[j-1 j+1])*G;
    end 
    Bi(j:j+1, j:m) = G' * Bi(j:j+1, j:m);
end
for j = 1:m-1
    Bi(1:j+1, j:j+1) = Bi(1:j+1, j:j+1) * G;
end
Bi(1,2) = 0;
end

%Calculate Givens rotation 2 x 2
function G = Givens_rotation( x ) 
    %Compute Givens rotation G so that G' * x = || x ||_2 % e_0        
    [ m, n ] = size( x );    
    assert( m==2 && n==1, 'x must be 2 x 1' );
    normx = norm( x );        
    gamma = x(1) / normx;
    sigma = x(2) / normx;  
    h = NaN;
    if(isequaln(gamma, h) >= 1)
        gamma = 0.00000000001;
    end
    if(isequaln(sigma, h) >=1 )
        sigma = 0.00000000001;
    end
    G = [ (gamma) (-sigma)
          (sigma)  (gamma) ]; 
end

