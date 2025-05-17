function [A,B,s] = truncqr(A,B,varargin)
    p = inputParser;
    addRequired(p,'A',@(x) ismatrix(x))
    addRequired(p,'B',@(x) ismatrix(x) && size(B,2) == size(A,2));
    addParameter(p,'tol',eps * size(A,2),@(x) isnumeric(x) && isscalar(x));
    parse(p,A,B,varargin{:});
    
    [Q,R] = qr(B,0);
    
    [U,S,V] = svd(A*R','econ');
    s = diag(S);
    rS = nnz(sqrt(cumsum(s.^2,'reverse')) > p.Results.tol * norm(s));
    U = U(:,1:rS);
    V = V(:,1:rS);
    S = S(1:rS,1:rS);
    
    A = U;
    B = Q*V*S;
end

