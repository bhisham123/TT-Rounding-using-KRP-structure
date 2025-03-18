function X = TTrounding_Two_Sided_Randomization_KRP(Y, l)

    [N,I,~] = TTsizes(Y);
    
    % Heuristic choice of oversampling ranks
    rho = [1; floor(2*l(2:N)'); 1];
    
    
    WL = KRPpartialContractionsLR(Y,max(l));
    WR = KRPpartialContractionsRL(Y, max(rho));
    
    
    X = cell(N,1);
    [U,S,V] = svd(WL{1}(1:l(2),:) * WR{1}(:,1:rho(2)), 'econ');
    S = diag(diag(S).^(-0.5));
    L = WR{1}(:,1:rho(2))*V*S;
    R = S*U'*WL{1}(1:l(2),:);
    
    X{1} = Y{1} * L;
    
    for n=2:N-1
       [U,S,V] = svd(WL{n}(1:l(n+1),:) * WR{n}(:,1:rho(n+1)), 'econ');
       S = diag(diag(S).^(-0.5));
       L = WR{n}(:,1:rho(n+1))*V*S;
       X{n} = h2v(R * v2h(Y{n} * L, I(n)), I(n));
       
       R = S*U'*WL{n}(1:l(n+1),:);
    end
    
    X{N} = h2v( R * v2h(Y{N}, I(N)), I(N));
end