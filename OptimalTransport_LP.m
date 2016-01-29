% solve optimal transport with linear programming
function [rcode, V, Pi, F] = OptimalTransport_LP(C, W1, W2)
% Solving the optimal transport problem in Euclidean space 
% Jianbo Ye <jxy198@ist.psu.edu> 2015
%
% - I/O
%   P1, P2 : n x d or m x d points
%   W1, W2 : n or m dimensional weights
%
% - Dependencies~
%   Mosek Matlab Toolbox: the network LP solver

n = size(W1,1);
m = size(W2,1);

global problemMap problemSet
key = [num2str(n)  'x'  num2str(m)];

if isKey(problemMap, key)
    prob = problemSet{problemMap(key)};
else
    subi = [repmat(1:n, 1, m), reshape(repmat((1:m) + n, n, 1), 1, n*m)];
    subj = [1:(n*m), 1:(n*m)];
    valij = ones(1,2*n*m);    
    prob.a = sparse(subi, subj, valij);    
    prob.blx = sparse(n*m, 1);
    prob.bux = [];
    problemSet{end+1} = prob;
    problemMap(key) = length(problemSet);
end

W1 = W1(:) / sum(W1(:));
W2 = W2(:) / sum(W2(:));

%assert(length(W1) == n && length(W2) == m);

% --- Mosek solver ---
prob.c = C; prob.c = prob.c(:);
prob.blc = [W1; W2];
prob.buc = prob.blc;

param.MSK_IPAR_OPTIMIZER = 'MSK_OPTIMIZER_PRIMAL_SIMPLEX'; 
[rcode, X] = mosekopt('minimize echo(0)', prob, param);


try
    V = dot(prob.c, X.sol.bas.xx);
    Pi = X.sol.bas.xx;
    F  = X.sol.bas.y;
    Pi(Pi<1E-10) = 0;
    Pi = sparse(reshape(Pi, n, m));
catch
    fprintf('MSKERROR: Could not get solution');
end



end