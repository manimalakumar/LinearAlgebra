function [Bi_next] = Bidiag_Francis_Step(Bi)
%  Input truly bi diagonal square matrix
%10.3.5.1
[m, m] = size(Bi);
%G1 w T's elements 
t00 = Bi(1,1) * Bi(1,1);
t10 = Bi(1,2) * Bi(1,1);
tm1m1 = ( Bi(m-1,m) * Bi(m-1,m) ) * ( Bi(m,m) * Bi(m,m));
G = eye(2);
IG = eye (m);
for j = 1:m-1
    if j == 1
        G = Givens_rotation( [ t00 - tm1m1;t10 ]);          
    else
        G= Givens_rotation([Bi(j,j);Bi(j+1,j)]);
    end 
    Bi(j:j+1, j:m) = G' * Bi(j:j+1, j:m);
end
for j = 1:m-1
    Bi(1:j+1, j:j+1) = Bi(1:j+1, j:j+1) * G;
end
Bi(1,2) = 0;
Bi_next = Bi;
end
    
%Calculate givens rotation 2 x 2
function G = Givens_rotation( x ) 
    %Compute Givens rotation G so that G' * x = || x ||_2 % e_0        
    [ m, n ] = size( x );    
    assert( m==2 && n==1, 'x must be 2 x 1' );
    normx = norm( x );        
    gamma = x(1) / normx;
    sigma = x(2) / normx;  
    
    G = [ (gamma) (-sigma)
          (sigma)  (gamma) ]; 
end

