function [ A_out, t_out, r_out ] = BiRed2( A, t, r )

  [ ATL, ATR, ...
    ABL, ABR ] = FLA_Part_2x2( A, ...
                               0, 0, 'FLA_TL' );

  [ tT, ...
    tB ] = FLA_Part_2x1( t, ...
                         0, 'FLA_TOP' );

  [ rT, ...
    rB ] = FLA_Part_2x1( r, ...
                         0, 'FLA_TOP' );
                     
  while ( size( ATL, 1 ) < size( A, 1 ) )

    [ A00,  a01,     A02,  ...
      a10t, alpha11, a12t, ...
      A20,  a21,     A22 ] = FLA_Repart_2x2_to_3x3( ATL, ATR, ...
                                                    ABL, ABR, ...
                                                    1, 1, 'FLA_BR' );

    [ t0, ...
      tau1, ...
      t2 ] = FLA_Repart_2x1_to_3x1( tT, ...
                                    tB, ...
                                    1, 'FLA_BOTTOM' );

    [ r0, ...
      rho1, ...
      r2 ] = FLA_Repart_2x1_to_3x1( rT, ...
                                    rB, ...
                                    1, 'FLA_BOTTOM' );
                                
    %------------------------------------------------------------%
    n = size(A);
    % Per HQR.m 3.3.4.2
    f = size (A22,1);
    if f >= 1 
        [rho, u21, tauL] = Housev(alpha11, a21);
        alpha11 = rho;
        a21 = u21;
        tau1 = tauL;
        w12t = (a12t + u21' * A22)/tauL; 
        a12t = a12t - w12t/tauL; % divide by tauL from UT paper
        A22 = A22 - u21 * w12t;          
    end
    % update row. Stop H(row) before H(col) Per 10.3.1.1       
    if f > 2 %A22 when not decomposed to 2 by 2 or 1 by 1
        [u12, tauR] = Housev1( a12t');
        a12t = u12';
        rho1 = tauR;
        u12(1, 1) = 1;           
        y21t = u21'*tril(A22,-1); % both sides H(u21,tau1) A22 H(u12,ρ1)
        B22 = A22 - (1/tauL)*u21*y21t;
        x12 = triu(B22,+1) * u12;
        A22 = B22 - (1/tauR)*x12*u12';       
    end

    if f < 1 %run only once when last row, end of while
        B = A;
        V = []; % lower Triangle w more H entries
        U = []; % upper Triangle w less H entries
        for k=1:n-1
         x = B(k:n,k);     
            % Householder reflection unit vector u from the vector x.
            m = max(abs(x));
            u = x/m;
            if u(1) == 0
             su = 1;
            else
             su = sign(u(1));
            end
            u(1) = u(1)+su*norm(u);
            u = u/norm(u);
            u = u(:);
         v = u;
         pos = k:n;
         B(pos,pos) = B(pos,pos) - 2*v*(v'*B(pos,pos)); %left householder reflection n projection
         v = [zeros(k-1,1);v];       
         V = [V v];     
         if k < n-1
             x = B(k,k+1:n)';          
                % Householder reflection unit vector u from the vector x.
                m = max(abs(x));
                u = x/m;
                if u(1) == 0
                 su = 1;
                else
                 su = sign(u(1));
                end
                u(1) = u(1)+su*norm(u);
                u = u/norm(u);
                u = u(:);
             row = 1:n;
             col = k+1:n;
             B(row,col) = B(row,col) - 2*(B(row,col)*u)*u'; %right householder orthogonal projection
             u = [zeros(k,1);u];
             U = [U u];            
         end
        end
    %Build the FLAME matrix Bidiagonal part. Keep orig householders, t/V, r/U.
    %Bidiagonal entries from new B
    A00 = B(1:n-1,1:n-1);
    n = n(1);
    a01 = B(1:n-1,n); 
    a10t = B(n,1:n-1); 
    alpha11 = B(n,n);   
    end 
    %------------------------------------------------------------%

    [ ATL, ATR, ...
      ABL, ABR ] = FLA_Cont_with_3x3_to_2x2( A00,  a01,     A02,  ...
                                             a10t, alpha11, a12t, ...
                                             A20,  a21,     A22, ...
                                             'FLA_TL' );

    [ tT, ...
      tB ] = FLA_Cont_with_3x1_to_2x1( t0, ...
                                       tau1, ...
                                       t2, ...
                                       'FLA_TOP' );

    [ rT, ...
      rB ] = FLA_Cont_with_3x1_to_2x1( r0, ...
                                       rho1, ...
                                       r2, ...
                                       'FLA_TOP' );
                                   
  end

  A_out = [ ATL, ATR
            ABL, ABR ];

  t_out = [ tT
            tB ];
        
  r_out = [ rT
            rB ];

return