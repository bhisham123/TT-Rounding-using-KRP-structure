function W = KRPpartialContractionsRL_Cookies(Y,X,n)
[N,I,~] = TTsizes(Y);
W = cell(N-1,1);
for i = N:-1:n+1
   if i == N
       W{i-1} = v2h(Y{i},I(i))*X{i-1};
   else
       W{i-1} = v2h(Y{i},I(i))*khatrirao(W{i},X{i-1});
   end
end
end