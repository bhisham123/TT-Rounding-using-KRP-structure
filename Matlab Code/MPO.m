% Heisenberg 1/2-spin system
% A = \sum_{i=1}^d X_{i} X_{i+1} + Y_{i} Y_{i+1} + Z_{i} Z_{i+1}

clear all; clc;
addpath("./TTcore/");
addpath("./TTrandomized/");
addpath("./cookies/");

d = 30;
J = 2^(-d/2);
X = sqrt(J)*[0,1;1,0]; X = X(:);
Y = sqrt(J)*[0,-1i;1i,0]; Y = Y(:);
Z = sqrt(J)*[1,0;0,-1]; Z = Z(:);
I = eye(2); I = I(:);

s = 3 * d - 3;
A = cell(s,1);
for k = 1 : 3
    if k == 1
        M = X;
    elseif k == 2
        M = Y;
    else
        M = Z;
    end

    for i = 1 : d - 1
        A{(k-1)*(d-1) + i} = cell(d,1);
        for j = 1 : d
            if j == i || j == i+1
                A{(k-1)*(d-1) + i}{j} = M;
            else
                A{(k-1)*(d-1) + i}{j} = I;
            end
        end
    end
end

t = tic;
T = A{1};
for i = 2 : s
    T = TTaxby(1,T,1,A{i});
end
fprintf("Actual norm: %e\n", TTnorm(T,"OLR"));

tol = 1e-8;
B = TTrounding(T,tol);
t = toc(t);
fprintf("Time for TT sum + TT SVD: %e\n", t);

t = tic;
B = TTsum_Randomize_then_Orthogonalize(A,ones(s,1),tol);
t = toc(t);
fprintf("Time for randTT sum: %e\n", t);


t = tic;
B = TTsum_Randomize_then_Orthogonalize_KRP(A, ones(s,1), tol);
t = toc(t);
fprintf("Time for randKRPTT sum: %e\n", t);


%%
% Operator from Camano, Epperly and Tropp preprint
clear all; clc;
d = 10; % in the paper this is chosen 101

s = d * d - d;

A = cell(s,1);
a = zeros(s,1);
alpha = 1.5;

J = 2^(-d/2);
X = sqrt(J)*[0,1;1,0]; X = X(:);
Y = sqrt(J)*[0,-1i;1i,0]; Y = Y(:);
Z = sqrt(J)*[1,0;0,-1]; Z = Z(:);
I = eye(2); I = I(:);

counter = 0;
for k = 1 : 2
    
    if k == 1
        M = X;
    else
        M = Z;
    end
    
    for i = 1 : d-1
        for j = i+1 : d
            counter = counter + 1;
            a(counter) = 0.5 * J/((j-i)^alpha);
            A{counter} = cell(d,1);
            for p = 1 : d
                if p == i || p == j
                    A{counter}{p} = M;
                else
                    A{counter}{p} = I;
                end
            end
        end
    end
end

t = tic;
T = A{1};
for i = 2 : s
    T = TTaxby(1,T,1,A{i});
end
T_norm = TTnorm(T);
fprintf("Actual norm: %e\n", T_norm);
tol = 1e-5;
B = TTrounding(T,tol);
t = toc(t);
fprintf("Time for TT sum + TT SVD: %e\n", t);
fprintf("Relative error-TT SVD: %e\n", TTnorm(TTaxby(1,T,-1,B))/T_norm);



t = tic;
B2 = TTsum_Randomize_then_Orthogonalize(A,ones(s,1),tol);
t = toc(t);
fprintf("Time for randTT sum: %e\n", t);
fprintf("Relative for randTT sum:: %e\n", TTnorm(TTaxby(1,T,-1,B2))/T_norm);

t = tic;
B3 = TTsum_Randomize_then_Orthogonalize_KRP(A, ones(s,1), tol/10);
t = toc(t);
fprintf("Time for randKRPTT sum: %e\n", t);
fprintf("Relative for randKRPTT sum:: %e\n", TTnorm(TTaxby(1,T,-1,B3))/T_norm);

