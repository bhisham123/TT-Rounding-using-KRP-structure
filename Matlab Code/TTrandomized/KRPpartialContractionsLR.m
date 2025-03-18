function W = KRPpartialContractionsLR(Y,rank)
[N,I,~] = TTsizes(Y);

% W = KTRrand(I, rank);
% X = cell(N-1,1);
% for i = 1:1:N-1
%    if i == 1
%        X{i} = W{i}'*Y{i};
%        % X{i-1} = v2h(Y{i},I(i))*W{i-1};
%        % Z{i-1} = X;
%    else
%        X{i} = khatrirao(W{i},X{i-1}')'*Y{i};
%        % X{i-1} = v2h(Y{i},I(i))*khatrirao(X{i},W{i-1});
%        % Z{i-1} = X;
%    end
% end

X = KTRrand(I, rank);
W = cell(N-1,1);
for i = 1:1:N-1
   if i == 1
       W{i} = X{i}'*Y{i};
   else
       W{i} = khatrirao(X{i},W{i-1}')'*Y{i};
   end
end

% W{1} = Y{1}'*X{1};
% for i = 2:1:N-1
%     W{i} = Y{i}'* khatrirao(X{i},W{i-1});
% end
% for i = 1:N-1
%     W{i} = W{i}';
% end


end