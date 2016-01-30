function [V, Pi] = OptimalTransport_IBP_Sinkhorn(C, W1, W2, rho, max_iter)

n = size(W1,1);
m = size(W2,1);

v=ones(m,1);

xi=exp(-C / rho);
%umax=zeros(max_iter);
%vmax=zeros(max_iter);
for i=1:max_iter
    u=W1 ./ (xi*v);
    v=W2 ./ (xi'*u);
%    umax(i)=max(u);
%    vmax(i)=max(v);
end

%figure; plot(log10(umax), log(vmax), 'd');
if any(isnan(u)) || any(isnan(v))
    error('NaN encountered.')
end
Pi=diag(u) * xi * diag(v);
V=sum(C(:).*Pi(:));