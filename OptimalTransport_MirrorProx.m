function [V, Pi] = OptimalTransport_MirrorProx(C, W1, W2, rho, max_iter)

n = size(W1,1);
m = size(W2,1);

assert(n==size(C,1) && m==size(C,2));

W2=W2';
Pi2=W1*W2;Pi1=Pi2;
Lambda=zeros(size(Pi2));
C=C/rho;
xi=exp(-C);

for i=1:max_iter    
    Lambda_c=Lambda;
    Pi1_c=Pi1;
    Pi2_c=Pi2;
    tmp = exp(Lambda);
    Lambda=Lambda_c + Pi1 - Pi2;
    Pi1=Pi2_c .* xi ./tmp  + eps;
    Pi1=bsxfun(@times, Pi1, W1 ./sum(Pi1, 2));
    Pi2=Pi1_c .* tmp + eps;
    Pi2=bsxfun(@times, Pi2, W2 ./sum(Pi2, 1));
    
    % false/true, mirror descent or mirror prox
    if (false)          
    tmp = exp(Lambda);
    Lambda=Lambda_c + Pi1 - Pi2;
    Pi1=Pi2_c .* xi ./tmp  + eps;
    Pi1=bsxfun(@times, Pi1, W1 ./sum(Pi1, 2));
    Pi2=Pi1_c .* tmp + eps;
    Pi2=bsxfun(@times, Pi2, W2 ./sum(Pi2, 1));    
    end
end

Pi=(Pi1 + Pi2)/2.;
V=rho*sum(C(:) .* Pi(:));
