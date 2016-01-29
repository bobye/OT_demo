function [V, Pi] = OptimalTransport_ABP_Sinkhorn(C, W1, W2, rho, max_iter)

n = size(W1,1);
m = size(W2,1);

v=ones(m,1);

xi=exp(-C / rho);
for i=1:max_iter
    u=W1 ./ (xi*v);
    v=W2 ./ (xi'*u);
end

Pi=diag(u) * xi * diag(v);
V=sum(C(:).*Pi(:));