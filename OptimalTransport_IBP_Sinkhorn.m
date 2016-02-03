function [V, Pi] = OptimalTransport_IBP_Sinkhorn(C, W1, W2, rho, max_iter)

n = size(W1,1);
m = size(W2,1);

W1=max(W1,eps); W1=W1/sum(W1);
W2=max(W2,eps); W2=W2/sum(W2);

%addpath('../AdvanpixMCT-3.9.1.10015');
%mp.Digits(34);
%v=mp(ones(m,1)); % quadruple precision
v=ones(m,1);

xi=exp(-C / rho);
if any(isnan(xi))
    error('NaN encountered.')
end
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
Pi=full(spdiags(u, 0, n, n) * xi * spdiags(v, 0, m, m));
Pi=double(Pi);
V=sum(C(:).*Pi(:));
V=double(V);