function [U,B,V] = Bidiag_Francis_Step_Update_U_V(U,B,V)
    % Iteration 1
    [m,~] = size(B);
    T = transpose(B)*B;
    x = [T(1,1)-T(m,m),T(2,1)];
    x = transpose(x);
    G = Givens_rotation(x);
    GI = eye(m,m);
    GI(1:2,1:2) = G;
    B = B * GI;
    V = V * GI;
    % Iteration 1 - 2nd part
    x = [B(1,1),B(2,1)];
    x = transpose(x);
    G = Givens_rotation(x);
    G = transpose(G);
    GI = eye(m,m);
    GI(1:2,1:2) = G;
    B = GI * B;
    U = U * GI;
    % Iteration 2
    for i=1:m-2
        x = [B(i,i+1),B(i,i+2)];
        x = transpose(x);
        G = Givens_rotation(x);
        GI = eye(m,m);
        GI(i+1:i+2,i+1:i+2) = G;
        B = B * GI;
        V = V * GI;
    % Iteration 2 - 2nd part
        x = [B(i+1,i+1),B(i+1,i+1)];
        x = transpose(x);
        G = Givens_rotation(x);
        G = transpose(G);
        GI = eye(m,m);
        GI(i+1:i+2,i+1:i+2) = G;
        B = GI * B;
        U = U * GI;
    end

end