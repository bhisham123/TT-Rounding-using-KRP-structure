function W = KRPpartialContractionsRL_New(Y,rank,n)
[N,I,~] = TTsizes(Y);

X = KTRrand_New(I, rank,n);
W = cell(N-1,1);
for i = N:-1:n+1
   if i == N
       W{i-1} = v2h(Y{i},I(i))*X{i-1};
       % Z{i-1} = X;
   else
       W{i-1} = v2h(Y{i},I(i))*khatrirao(W{i},X{i-1});
       % Z{i-1} = X;
   end
end
end