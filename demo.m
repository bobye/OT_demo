close all
%% transport between two Gaussian Mixtures
N=128;
p1 = gmdistribution([0.2; 0.5], reshape([0.05, 0.03].^2, 1, 1, 2), [0.5, 0.5]);
p2 = gmdistribution([0.6; 0.8], reshape([0.05, 0.05].^2, 1, 1, 2), [0.4, 0.6]);

locations=[0:1/(N-1):1]';
w1=p1.pdf(locations) + 0.1; w1=w1(:)/sum(w1(:));
w2=p2.pdf(locations) + 0.1; w2=w2(:)/sum(w2(:));

C=pdist2(locations, locations).^2;
plot(locations, w1); hold on
plot(locations, w2);


max_iters=[1, 10, 50, 200, 1000, 5000];


%% 
OptimalTransport_LPServer;
figure;

for s=1:6
%% linear programming
[~, V1, Pi1, ~] = OptimalTransport_LP(C, w1, w2);
Pi1=full(Pi1); 
%subplot_tight(4,6,1+(s-1)*5); imshow(-Pi1, []);


%% iterative Bregman projection: Sinkhorn distance
tic;
[V2, Pi2] = OptimalTransport_IBP_Sinkhorn(C, w1, w2, .1/N, max_iters(s));
toc;
subplot_tight(6,4,1+(s-1)*4); imshow(addframe(imfuse(-Pi2, -Pi1)));

[V2, Pi2] = OptimalTransport_IBP_Sinkhorn(C, w1, w2, .5/N, max_iters(s));
subplot_tight(6,4,2+(s-1)*4); imshow(addframe(imfuse(-Pi2, -Pi1)));

[V2, Pi2] = OptimalTransport_IBP_Sinkhorn(C, w1, w2, 2./N, max_iters(s));
subplot_tight(6,4,3+(s-1)*4); imshow(addframe(imfuse(-Pi2, -Pi1)));

%% Bregman ADMM
rho=2.*mean(C(:));
tic;
[V3, Pi3] = OptimalTransport_BADMM(C, w1, w2, rho, max_iters(s));
toc;
subplot_tight(6,4,4+(s-1)*4); imshow(addframe(imfuse(-Pi3, -Pi1)));


end
tightfig;


