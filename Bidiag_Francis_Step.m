function [BiNext] = Bidiag_Francis_Step( Bi )
    % Iteration 1
    [m,~] = size(Bi);
    T = transpose(Bi)*Bi;
    x = [T(1,1)-T(m,m),T(2,1)];
    x = transpose(x);
    G = Givens_rotation(x);
    GI = eye(m,m);
    GI(1:2,1:2) = G;
    Bi = Bi * GI;
    % Iteration 1 - 2nd part
    x = [Bi(1,1),Bi(2,1)];
    x = transpose(x);
    G = Givens_rotation(x);
    G = transpose(G);
    GI = eye(m,m);
    GI(1:2,1:2) = G;
    Bi = GI * Bi;
    % Iteration 2
    for i=1:m-2
        x = [Bi(i,i+1),Bi(i,i+2)];
        x = transpose(x);
        G = Givens_rotation(x);
        GI = eye(m,m);
        GI(i+1:i+2,i+1:i+2) = G;
        Bi = Bi * GI;
    % Iteration 2 - 2nd part
        x = [Bi(i+1,i+1),Bi(i+1,i+1)];
        x = transpose(x);
        G = Givens_rotation(x);
        G = transpose(G);
        GI = eye(m,m);
        GI(i+1:i+2,i+1:i+2) = G;
        Bi = GI * Bi;
    end
    BiNext = Bi;
end
