function Y = LRLGramRoundTT(X, tau, varargin)
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
    
    
    Y = X;
    GR = cellfun(@(x) (x + x')/2, nestedContractionsRL(X,X), 'UniformOutput', false);
    
    tol = tau / sqrt(N-1);
    % tol = tau * sqrt(GR{1}) / sqrt(N - 1); 
    
    temp = Y{1}' * Y{1};
    temp = (temp + temp')/2; 
    [ML, MR] = truncgram(temp, GR{2},'tol',tol,'alg',p.Results.alg);
    Y{1} = Y{1} * ML;
    Y{2} = H2V(MR' * V2H(Y{2}, I(2)), I(2)); 
    for i = 2 : N - 1
        temp = Y{i}' * Y{i};
        temp = (temp + temp')/2; 
        [ML, MR] = truncgram(temp, GR{i + 1},'tol',tol,'alg',p.Results.alg);
        Y{i} = Y{i} * ML;
        Y{i + 1} = H2V(MR' * V2H(Y{i + 1}, I(i + 1)), I(i + 1));
    end
end