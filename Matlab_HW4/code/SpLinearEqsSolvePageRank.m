function ranking = SpLinearEqsSolvePageRank(G)

N = size(G, 1);
c = sum(G, 1);
[I,J,K] = find(c);
K = 1./K;
d = sparse(I, J, K);
D = diag(d);

rho = 0.85;
e = ones(N,1);
I = speye(N);
PR = (I-rho*G*D)\e;

ranking = sortrows([PR, (1:N)'], -1);
ranking = ranking(:,2);
end