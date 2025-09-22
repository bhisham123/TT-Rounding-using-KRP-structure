function X = TTsum_Randomize_then_Orthogonalize_KRP(S, a, tol, incr_f, samples)
%TTsum_Randomize_then_Orthogonalize_KRP: Randomized TT-sum with Khatri-Rao projections and adaptive rounding.
%
%   X = TTsum_Randomize_then_Orthogonalize_KRP(S, a, tol, incr_f, use_svd_filtering, svd_thresh)
%   computes the weighted sum of a set of tensors in Tensor-Train (TT) format using
%   randomized projections based on the Khatri-Rao product and an adaptive rank expansion
%   strategy. 
%
%   Inputs:
%     - S                : Cell array of TT-tensors to be summed. Each element is a TT-format tensor.
%     - a                : Vector of weights corresponding to each tensor in S (same length as S).
%     - tol              : Desired relative approximation tolerance (Frobenius norm). Default: 1e-10.
%     - incr_f           : Increment factor controlling the block size for randomized projections.
%                          The block size is computed as floor(maxr*incr_f). Default: 0.05.
%     - smaples          : minimum number of smaples to use for norm estimation, Default: 20.
%     
%   Output:
%     - X : The TT-format tensor (as a cell array of cores) representing the approximate weighted sum
%           of the input TT-tensors in S.
%

% ----------- Default Parameters --------------
if nargin < 5 || isempty(samples)
    samples = 20;
end
if nargin < 4 || isempty(incr_f)
    incr_f = 0.05;
end
if (nargin == 2)
    tol = 1e-10;
end


m = length(S);
assert(length(a)==m, "Incompatible size of coefficients and summands in TTrounding_Sum");

[N,I,r1] = TTsizes(S{1});
rs = zeros(N+1,m);
rs(:,1) = r1;
for j=2:m
   rs(:,j) = TTranks(S{j}); 
end

max_mod_rank = max(rs, [], 2);

for j = 1:m
    S{j} = TTscale(S{j}, a(j));
end


% Form the complete random projections
W = cell(N-1,1);


temp_init_b = max(max(max_mod_rank),samples); %initial block size
O = KRPrand(I, temp_init_b, 1);
W_Ind_New = cell(N-1,m);
for j = 1:m
    W_Ind_New(:,j) = KRPpartialContractionsRL_Cookies(S{j}, O, 1);
end
for j = 1:N-1
    W{j} = [W{j},vertcat(W_Ind_New{j,:})];
end


% Initialize with the first core of the formal TT-sum
X = cell(N, 1);
lr = [0, cumsum(rs(2, :))];
X{1} = zeros(I(1), lr(m+1));
for j = 1:m
    X{1}(:, lr(j)+1 : lr(j+1)) = S{j}{1};
end

normX_est = norm(X{1}*W{1},'fro')/sqrt(temp_init_b); %Estimate norm of the sum
tau = tol*normX_est/sqrt(N-1);  %error distributed to each core

% Left-to-right randomization and orthogonalization
for n = 1:N-1
    current_core = X{n}; %current core of the TT-sum
    [rows,cols] = size(current_core);
    cols = min(rows,cols);
    Init_b = max_mod_rank(n+1); %initial block size
    orthogonal_cols = 0;

    %compute partial sketches
    if size(W{n},2) < Init_b
        O = KRPrand(I, Init_b-size(W{n},2), n);
        W_Ind_New = cell(N-1,m);
        for j = 1:m
            W_Ind_New(:,j) = KRPpartialContractionsRL_Cookies(S{j}, O, n);
        end
        for j = n:N-1
            W{j} = [W{j},vertcat(W_Ind_New{j,:})];
        end
    end


    Yn = current_core * W{n}(:,1:Init_b); %compute sketch

    [X{n}, ~] = qr(Yn, 0); %tall and skinny QR

    % Pass the remaning factor to core on right
    Mn = X{n}' * current_core;
    lr = [0, cumsum(rs(n + 1, :))];
    rr = [0, cumsum(rs(n + 2, :))];
    if(n < N - 1)
        X{n+1} = zeros(size(Mn,1) * I(n + 1), rr(m + 1));
        for j = 1:m
            x =  Mn(:, lr(j)+1:lr(j+1)) * v2h(S{j}{n+1}, I(n+1));
            X{n+1}(:, rr(j)+1:rr(j+1)) = h2v(x, I(n+1));
        end
    else
        X{N} = Mn(:, 1:lr(2)) * v2h(S{1}{N}, I(N));
        for j = 2:m
            X{N} = X{N} + Mn(:, lr(j)+1:lr(j+1)) * v2h(S{j}{N}, I(N));
        end
        X{N} = h2v(X{N}, I(N));
    end
    

    orthogonal_cols = orthogonal_cols + Init_b; %update the orthogonal columns count    

    b = max(1,floor(cols*incr_f)); %block size for increment

    tempb = max(b,samples); %samples for norm estimation

    %compute partial sketches
    if size(W{n},2) < orthogonal_cols+tempb
        O = KRPrand(I, (orthogonal_cols+tempb)-size(W{n},2), n);
        W_Ind_New = cell(N-1,m);
        for j = 1:m
            W_Ind_New(:,j) = KRPpartialContractionsRL_Cookies(S{j}, O, n);
        end
        for j = n:N-1
            W{j} = [W{j},vertcat(W_Ind_New{j,:})]; 
        end
    end
    Yn = current_core*W{n}(:,orthogonal_cols+1:orthogonal_cols+tempb); %compute sketch
    Yn = Yn - X{n}*(X{n}'*Yn); % Remove the directions already computed (orthogonalization)

    while (norm(Yn,'fro')/sqrt(tempb) > tau) 
        b = min(b,cols - orthogonal_cols); 
        if b <= 0
            break
        end

        Yn = Yn(:,1:b);
        [Q, ~] = qr(Yn, 0);   %tall and skinny QR
        [Q, ~] = qr(Q-X{n}*(X{n}'*Q), 0); % Reorthogonalization
       
        X{n} = [X{n},Q];  % Add new directions

        %Pass remaining factor to the core on right
        Mn =  Q'*current_core;
        if(n < N - 1)
            TempX = zeros(size(X{n},2) * I(n + 1), rr(m + 1));
            for j = 1:m
                x =  Mn(:, lr(j)+1:lr(j+1)) * v2h(S{j}{n+1}, I(n+1));
                TempX(:, rr(j)+1:rr(j+1)) = h2v([v2h(X{n+1}(:, rr(j)+1:rr(j+1)),I(n+1));x], I(n+1));
            end
            X{n+1} = TempX;
        else
            TempX = Mn(:, 1:lr(2)) * v2h(S{1}{N}, I(N));
            for j = 2:m
                TempX = TempX + Mn(:, lr(j)+1:lr(j+1)) * v2h(S{j}{N}, I(N));
            end
            X{N} = h2v([v2h(X{N},I(N));TempX], I(N));
        end

        orthogonal_cols = orthogonal_cols + b;

        tempb = max(b,samples);

        %Compute partial skeches
        if size(W{n},2) < orthogonal_cols+tempb
            O = KRPrand(I, (orthogonal_cols+tempb)-size(W{n},2), n);
            W_Ind_New = cell(N-1,m);
            for j = 1:m
                W_Ind_New(:,j) = KRPpartialContractionsRL_Cookies(S{j}, O, n);
            end
            for j = n:N-1
                W{j} = [W{j},vertcat(W_Ind_New{j,:})];
            end
        end

        Yn = current_core*W{n}(:,orthogonal_cols+1:orthogonal_cols+tempb); %compute sketch
        Yn = Yn - X{n}*(X{n}'*Yn); % Remove the directions already computed (orthogonalization)
    end
    if norm(X{n}'*X{n} - eye(size(X{n},2)),'fro') > 1e-7
        disp([' n = ', num2str(n), ', Othogonality Error ', num2str(norm(X{n}'*X{n} - eye(size(X{n},2)),'fro')), ', Condition number = ', num2str(cond(current_core),'%.2e')])
        error('Orthonormality Error')
    end
end

% TT-rounding of X, which is already orthogonalized.
normX = norm(X{N}, 'fro');
tau = tol * normX / sqrt(N - 1);                    % Truncation threshold
for n = N:-1:2
    [Q,R] = qr(v2h(X{n}, I(n))', 0);
    [U, S, V] = svd(R, 'econ');
    rank = trunc(diag(S), tau);
    X{n} = h2v( (Q*U(:, 1:rank))', I(n));
    try
    X{n-1} = X{n-1} * (V(:, 1:rank) * S(1:rank, 1:rank));
    catch
        warning('Some problem');
        break;
    end
end
end


