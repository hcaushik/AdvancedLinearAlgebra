function [ A_out, t_out, r_out ] = BiRed( A, t, r )

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
    [m, ~] = size(ATL);
    [m2,~] = size(A);
    if (m > 0)
        x = [A00(m,m), a10t(1,m), transpose(A20(:,m))];
        x = transpose(x);
        [x,tau1] = Housev1(x);
        A00(m,m) = x(1);
        a10t(1,m) = x(2);
        A20(:,m) = x(3:end);
        Aright(1,1:m2-m) = [a01(m),A02(m,:)];
        Aright(2,1:m2-m) = [alpha11, a12t];
        Aright(3:m2-(m-1), 1:m2-m) = cat(2,a21(1:m2-m-1),A22);
        I = eye(m2-m+1);
        u21 = [1,transpose(x(2:end))];
        u21 = transpose(u21);
        house_transform = I - (1/tau1)*(u21)*transpose(u21);
        Aright = house_transform * Aright;
        if m < m2
            [u12, rho1] = Housev1(transpose(Aright(1,:)));
            Aright(1,:) = transpose(u12);
            u12(1,1) = 1;
            I = eye(m2-m);
            house_transform = I - (1/rho1)*u12*transpose(u12);
            Aright(2:end,:) = Aright(2:end,:) * house_transform;
            size([a01(m),A02(m,:)]);
            size(Aright(1,1:m2-m));
            a01(m) = Aright(1,1);
            A02(m,:) = Aright(1,2:end);
            alpha11 = Aright(2,1);
            a12t = Aright(2, 2:m2-m);
            a21 = Aright(3:m2-(m-1), 1);
            A22 = Aright(3:m2-(m-1), 2:m2-m);
            clear Aright;
        end
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