function [V, Pi] = OptimalTransport_BADMM(C, W1, W2, rho, max_iter)

n = size(W1,1);
m = size(W2,1);

Pi2=W1*W2';
Lambda=zeros(size(Pi2));
C=C/rho;
xi=exp(-C);
for i=1:max_iter
    tmp = exp(Lambda);
    Pi1=Pi2 .* xi ./tmp  + eps;
    Pi1=bsxfun(@times, Pi1, W1 ./sum(Pi1, 2));
    Pi2=Pi1 .* tmp + eps;
    Pi2=bsxfun(@times, Pi2', W2 ./sum(Pi2, 1)')'; % memory overheads in matrix transpose
    Lambda=Lambda + Pi1 - Pi2;
end

Pi=(Pi1 + Pi2)/2.;
V=rho*sum(C(:) .* Pi(:));
