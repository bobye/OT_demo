function [V, Pi, F] = OptimalTransport_Simulated_Annealing(C, W1, W2, rho, max_iter)

n = size(W1,1);
m = size(W2,1);
W1=max(W1,eps); W1=W1/sum(W1);
W2=max(W2,eps); W2=W2/sum(W2);

f1 = zeros(n,1);
T=rho;
obj=zeros(max_iter,1);
sigma=power(1/max_iter^3, 1/max_iter);
for i=1:max_iter
    L = max(bsxfun(@plus, -C, f1), [], 1);
    f2 = L + exprnd(T./W2' , [1, m]);
    U = min(bsxfun(@plus,  C, f2), [], 2);
    f1 = U - exprnd(T./W1, [n, 1]);
    T=T*sigma;    
    obj(i)=dot(f1, W1) - dot(f2, W2);
end
%close all
%plot(obj);%ylim([0, 0.2]);
F=[U(:); L(:)];  
V=dot(F, [W1; -W2]);

% method 1
if true
[~, midx2] = max(bsxfun(@plus,-C, U), [], 1);
[~, midx1] = min(bsxfun(@plus, C, L), [], 2);
Pi= zeros(n,m);
midx1=sub2ind([n,m], 1:n, midx1');
midx2=sub2ind([n,m], midx2, 1:m);
Pi(midx1)=W1; Pi(midx2)=Pi(midx2)+W2'; Pi=Pi/2;
end

% method 2
if false
s2=exp(-bsxfun(@plus, C, L)/T); s2=bsxfun(@times, s2, 1./sum(s2,2));
s1=exp(-bsxfun(@plus,-C, U)/T); s1=bsxfun(@times, s1, 1./sum(s1,1));
s2=bsxfun(@times, s2, W1);
s1=bsxfun(@times, s1, W2);
Pi=s1 + s2;
Pi=Pi/ sum(Pi(:));
end


%figure;
%plot(U(:)); hold on;
%plot(L(:)); hold on;
%plot(f1(:), 'o'); hold on;
%plot(f2(:), 'o');
end