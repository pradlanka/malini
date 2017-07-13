function Y=multinomial(x,W,N,P)

Y=zeros(N,P);

sum=ones(N,1);
for p=1:P
    size(exp(x{p}*W{p}))
    sum=sum+ exp(x{p}*W{p});
end

for p=1:P
    Y(:,p)=exp(x{p}*W{p})./sum;
end