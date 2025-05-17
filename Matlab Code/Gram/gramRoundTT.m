function Y = gramRoundTT(X, tau, varargin)
    p = inputParser;
    addRequired(p,'X',@(x) iscell(x) && isvector(X));
    addRequired(p,'tau',@(x) isscalar(x) && x >= 0);
    addParameter(p,'alg','eigen',@(x) ischar(x) && any(strcmp(x, {'eigen','chol'})));
    addParameter(p,'ordered',false,@(x) isscalar(x));
    parse(p,X,tau,varargin{:});
    
    H2V = @(X, I) reshape(X, [size(X, 1) * I, size(X, 2)/I]);
    V2H = @(X, I) reshape(X, [size(X, 1)/I, I * size(X, 2)]);
    
    N = numel(X);
    
    I = zeros(N, 1);
    I(1) = size(X{1}, 1);
    for i=2:N
        I(i) = size(X{i}, 1)/size(X{i - 1}, 2);
    end
    
    tol = tau / sqrt(N-1);
    
    GL = cellfun(@(x) (x + x')/2, nestedContractionsLR(X,X), 'UniformOutput', false);
    GR = cellfun(@(x) (x + x')/2, nestedContractionsRL(X,X), 'UniformOutput', false);
    
    ML = cell(N,1);
    MR = cell(N,1);
    
    [ML(1:N-1), MR(2:N)] = cellfun(@(A,B) truncgram(A,B,'tol',tol,'alg',p.Results.alg), GL(1:N-1), GR(2:N), 'UniformOutput', false);

    Y = X;
    Y{1} = X{1} * ML{1};
    for i=2:N-1
        Y{i} = H2V(MR{i}' * V2H(X{i} * ML{i},I(i)),I(i));
    end
    Y{N} = H2V(MR{N}' * V2H(X{N},I(N)),I(N));
end