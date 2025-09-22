function X = TTroundingKRP_Adaptive(Y,tol,init_f,incr_f,tol_scale,samples,spherical,do_rounding)
%TTroundingKRP_Adaptive: Adaptive TT-rounding using Khatri-Rao-based projections.
%
%   X = TTroundingKRP_Adaptive(Y, tol, init_f, incr_f, use_svd_filtering)
%   performs adaptive Tensor-Train (TT) rounding on a TT-format tensor Y 
%   using randomized projections based on the Khatri-Rao product and an 
%   incremental block expansion strategy.
%
%   Inputs:
%     - Y                 : A TT-tensor represented as a cell array of cores 
%                          (each vertically unfolded). Size N-by-1.
%     - tol               : Desired approximation accuracy in Frobenius norm.
%     - init_f (optional) : Initial factor to compute the starting block size 
%                          for projections, i.e., b_init = floor(R(n+1)*init_f).
%                          The default value of init_f is 0.1.
%     - incr_f (optional) : Incremental factor used to compute the block size of
%                          incremental QR, i.e., b = floor(R(n+1)*incr_f).
%                          The default value of incr_f is 0.05.
%     - samples (optional): minimum number of smaples to use for norm estimation, Default: 20.
%     - spherical (optional): Binary variable. True means use spherical distribution and false means Gaussian.
%                             Default value is false.
%     - do_rounding (optional):Binary variable. True means perform rounding pass at the end.
%                              Default value is false.
%

%
%   Output:
%     - X : A TT-tensor (cell array of cores) with reduced TT-ranks, 
%           approximating the original tensor Y up to the specified tolerance.
%
%   Description:
%     The function iteratively projects each TT core onto an orthonormal basis
%     constructed using randomized range finding. It begins with a small number
%     of random projections determined by `init_f`, and adaptively increases the 
%     projection dimension using `incr_f` until the approximation error falls 
%     below a locally computed tolerance.
%
%     The contraction matrices used in the Khatri-Rao product computations are 
%     updated in blocks using 'KRPpartialContractionsRL_New' when needed.

% ----------- Default Parameters --------------
if nargin < 8 || isempty(do_rounding)
    do_rounding = false; 
end
if nargin < 7 || isempty(spherical)
    spherical = false; 
end
if nargin < 6 || isempty(samples)
    samples = 20; 
end
if nargin < 5 || isempty(tol_scale)
    tol_scale = 1; 
end
if nargin < 4 || isempty(incr_f)
    incr_f = 0.05;
end
if nargin < 3 || isempty(init_f)
    init_f = 0.1;
end

[N,I,r] = TTsizes(Y);
W = cell(N-1,1);

X = Y;

% filename = [num2str(tol) '.txt'];
% fileID = fopen(filename, 'a');
% fprintf(fileID, 'Ranks:');
% fprintf(fileID, ' %2d', r);
% fprintf(fileID, '\n');

% Temp_Init_b = ones(N+1,1);
% Temp_Init_b(2:N) = floor(r(2:N)*init_f);
% Temp_Init_b(2) = max(Temp_Init_b(2),20);

temp_init_b = max(max(floor(r(2:N)*init_f)),samples);
W = KRPPartialContractionsRL(Y, temp_init_b,1,spherical);
normX_est = norm(X{1}*W{1},'fro')/sqrt(temp_init_b); % tensor norm estimation
tau = tol*normX_est/sqrt(N-1);  %error distributed to each core

% actual_norm = TTnorm(Y);
% fprintf(fileID, 'Actual norm of input tensor: %.10f', actual_norm);
% fprintf(fileID, '\tEstimate: %e', normX_est);
% fprintf(fileID, '\tRelative Error: %e', abs(normX_est-actual_norm)/actual_norm);
% fprintf(fileID, '\t Actual Tau: %e\n\n', tol*actual_norm/sqrt(N-1));

%Rouding tensor cores from left to right
for n = 1 : N - 1
    Current_core = X{n};
    [rows,cols] = size(Current_core);
    max_rank = min(rows,cols);
    Next_core    = X{n+1};
    
    Init_b = max(floor(max_rank*init_f),1); %initial block size
    orthogonal_cols = 0; 
    
   
    if size(W{n},2) < Init_b
         W_new = KRPPartialContractionsRL(Y,Init_b-size(W{n},2),n,spherical);
        for j = n:N-1
            W{j} = [W{j},W_new{j}];
        end
    end
    
    S = Current_core*W{n}(:,1:Init_b); %Compute sketch

    [X{n}, ~] = qr(S, 0); %tall and skinny QR 
    
    %compute the low rank approximation 
    M = X{n}' * Current_core;
    X{n+1} =  h2v(M * v2h(Next_core, I(n+1)), I(n+1));
    
    orthogonal_cols = orthogonal_cols + Init_b; %update the orthogonal columns count
    
    
    
    b = max(floor(max_rank*incr_f),1); %block size for increment
    tempb = max(b,samples);  %samples for norm estimation
    
    % compute partial sketches
    if size(W{n},2) < orthogonal_cols + tempb
        W_new = KRPPartialContractionsRL(Y,(orthogonal_cols + tempb)-size(W{n},2),n,spherical);
        for j = n:N-1
            W{j} = [W{j},W_new{j}];
        end
    end
    
    S = Current_core*W{n}(:,orthogonal_cols+1:orthogonal_cols+tempb); %compute sketch
    S = S - X{n}*(X{n}'*S); % Remove the directions already computed (orthogonalization)

    % fprintf(fileID, 'n = %2d, ', n);
    % actual_resd_norm = TTnorm(TTaxby(1,Y,-1,X));
    % fprintf(fileID, 'Actual residual norm: %.10f', actual_resd_norm);
    % resd_norm_est = norm(S,'fro')/sqrt(b);
    % fprintf(fileID, '\tEstimate: %.10f', resd_norm_est);
    % fprintf(fileID, '\tRelative Error: %.10f\n', abs(resd_norm_est - actual_resd_norm)/actual_resd_norm);
    
    while (norm(S,'fro')/sqrt(tempb) > tau/tol_scale) 
        b = min(b,max_rank - orthogonal_cols); 
        if b <= 0
            break
        end
        
        S = S(:,1:b);
        [Q,~] = qr(S, 0); %tall and skinny QR
        
        %Internal Reorthogonalization
        [Q,~] = qr(Q - X{n}*(X{n}'*Q),0);
    
        X{n} = [X{n},Q]; %Add new directions
    
        M = Q' * Current_core;
        V = M * v2h(Next_core, I(n+1));
        X{n+1} =  h2v([v2h(X{n+1}, I(n+1));V], I(n+1));
    
        orthogonal_cols = orthogonal_cols + b; %update the orthogonal columns count
    
        tempb = max(b,samples); %samples for norm estimation

        % Compute partial sketches
        if size(W{n},2) < orthogonal_cols+tempb
            W_new = KRPPartialContractionsRL(Y,(orthogonal_cols + tempb)-size(W{n},2),n,spherical);
            for j = n:N-1
                W{j} = [W{j},W_new{j}];
            end
        end
        S = Current_core*W{n}(:,orthogonal_cols+1:orthogonal_cols+tempb); %compute sketch
        S = S - X{n}*(X{n}'*S); % Remove the directions already computed (orthogonalization)

        % fprintf(fileID, 'n = %2d, ', n);
        % actual_resd_norm = TTnorm(TTaxby(1,Y,-1,X));
        % fprintf(fileID, 'Actual residual norm: %.10f', actual_resd_norm);
        % resd_norm_est = norm(S,'fro')/sqrt(b);
        % fprintf(fileID, '\tEstimate: %.10f', resd_norm_est);
        % fprintf(fileID, '\tRelative Error: %.10f\n', abs(resd_norm_est - actual_resd_norm)/actual_resd_norm);
    
    end
     if norm(X{n}'*X{n} - eye(size(X{n},2)),'fro') > 1e-7
        disp([' n = ', num2str(n), ', Othogonality Error ', num2str(norm(X{n}'*X{n} - eye(size(X{n},2)),'fro')), ', Condition number = ', num2str(cond(Current_core),'%.2e')])
        disp('')
     end
end

% [~,~,r1] = TTsizes(X);
% fprintf(fileID, 'Compressed ranks:');
% fprintf(fileID, ' %2d', r1);
% fprintf(fileID, '\n');
% 
% fprintf(fileID, 'Relative error:');
% fprintf(fileID, ' %.10f', TTnorm(TTaxby(1,X,-1,Y), "OLR") / actual_norm);
% fprintf(fileID, '\n\n');



% % TT-rounding of X, which is already orthogonalized.
if do_rounding
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
end



