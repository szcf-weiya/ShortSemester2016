function [iterations, ranking] = IterationSolvePageRank(G,error)
N = size(G, 1); % number of website
c = sum(G, 1);
ranking = [];
[~, index] = find(c > 0);
for i = index
    G(:, i) = G(:, i)./c(i);
end
iterations = 0;
PR0 = ones(N, 1);
PR = G * PR0;
while norm(PR - PR0,2) > error
	PR0 = PR;
	PR = G * PR0;
	iterations = iterations + 1;
    if iterations > 10000
        disp('`10000步迭代后无法达到指定误差限`');
        return;
    end
end
message = [num2str(iterations), '`步后达到指定误差限`'];
disp(message);
ranking = sortrows([PR,(1:N)'], -1);
ranking = ranking(:,2);
end