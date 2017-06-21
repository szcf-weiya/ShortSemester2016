function ranking = LinearEqsSolvePageRank(G)

N = size(G, 1);
c = sum(G, 1);
rho = 0.85;
d = zeros(N, 1);
[~,index] = find(c > 0);
d(index) = 1./c(index);

D = diag(d);
e = ones(N,1);
I = eye(N);
PR = (I - rho*G*D)\e;

ranking = sortrows([PR, (1:N)'], -1);
ranking = ranking(:,2);
end
