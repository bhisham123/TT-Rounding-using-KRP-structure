function [MA,MB,s] = truncgram(AtA,BtB,varargin)
    p = inputParser;
    addRequired(p,'AtA',@(x) ismatrix(x))
    addRequired(p,'BtB',@(x) ismatrix(x) && all(size(x) == size(AtA)));
    addParameter(p,'tol',eps * size(AtA,2),@(x) isnumeric(x) && isscalar(x));
    addParameter(p,'alg','eigen',@(x) ischar(x) && any(strcmp(x, {'eigen','chol'})));
    parse(p,AtA,BtB,varargin{:});
    
    switch p.Results.alg
        case 'eigen'
            [VA,DA] = eig(AtA, 'vector');
            [VB,DB] = eig(BtB, 'vector');
               
            % negate singular vectors based on negative eigenvalues
            VA(:,DA < 0) = -VA(:,DA < 0);
            VB(:,DB < 0) = -VB(:,DB < 0);
            
            % absolute value of eigenvalues and sort
            [SA, IA] = sort(abs(DA),'descend');
            [SB, IB] = sort(abs(DB),'descend');
            
            % reorder singular vectors
            VA = VA(:,IA);
            VB = VB(:,IB);
            
            % compute singular values from eigenvalues
            SA = diag(sqrt(SA));
            SB = diag(sqrt(SB));
            
            % A *= VA*pinv(SA)*U
            MA = VA*pinv(SA);
            % B *= VB*pinv(SB)*(V*S)
            MB = VB*pinv(SB);
            
            % M = SA*VA'*VB*SB;
            M = SA*VA'*VB*SB;
        case 'chol'
            [RA,pA] = cholpivot(AtA); % AtA(pA,pA) = RA'*RA
            [RB,pB] = cholpivot(BtB); % BtB(pB,pB) = RB'*RB
            
            % truncate zero rows to get R1 = [R11 R12] blocks
            RA1 = RA(any(RA,2),:);
            RB1 = RB(any(RB,2),:);

            % invert permutation vectors
            pA(pA) = 1:length(pA); 
            pB(pB) = 1:length(pB); 

            % A *= PA'*pinv(RA)*U
            %MA = pinv(RA(:,pA));
            % B *= PB'*pinv(RB)*(V*S)
            %MB = pinv(RB(:,pB));
            
            % M = RA*PA*PB'*RB'
            M = RA1(:,pA)*RB1(:,pB)';
    end
    
    [U,S,V] = svd(M,'econ');
    s = diag(S);
    rS = nnz(sqrt(cumsum(s.^2,'reverse')) > p.Results.tol * norm(s));
    U = U(:,1:rS);
    V = V(:,1:rS);
    S = S(1:rS,1:rS);
    
    switch p.Results.alg
        case 'eigen'
            MA = MA*U;
            MB = MB*(V*S);
        case 'chol'
            % truncate columns to get R11 block
            rnkA = size(RA1,1);
            RA11 = RA1(:,1:rnkA);
            rnkB = size(RB1,1);
            RB11 = RB1(:,1:rnkB);
            
            % MA = PA'*[inv(RA11); 0]*U
            MA = zeros(length(pA),size(U,2));
            MA(1:rnkA,:) = RA11\U;
            MA = MA(pA,:);
            
            % MB = PB'*[inv(RB11); 0]*V*S
            MB = zeros(length(pB),size(V,2));
            MB(1:rnkB,:) = RB11\(V*S);
            MB = MB(pB,:);
    end
    
    s = [s; NaN(size(AtA,1)-length(s),1)];
end