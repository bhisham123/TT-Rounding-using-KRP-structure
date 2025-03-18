function X = TTKrpSketch_old(Y,W,rank)
[N,I,~] = TTsizes(Y);
% X = zeros(size(Y{1},1),size(W{1},2));
for i = N:-1:2
   if i == N
       X = v2h(Y{i},I(i))*W{i-1}(:, 1:rank);
   else
       X = v2h(Y{i},I(i))*khatrirao(X,W{i-1}(:, 1:rank));
   end
end
X = Y{1}*X;
end