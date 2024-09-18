% This script runs a test that runs Monte-Carlo Markov Chains simulation
% for the problem

% Adjacency matrix for the graph
G = [
0 1 1 1 0 0 0 0;
1 0 1 0 0 1 0 0;
1 1 0 1 0 0 0 0;
1 0 1 0 0 0 0 1;
0 0 0 0 0 1 1 1;
0 1 0 0 1 0 1 0;
0 0 0 0 1 1 0 1;
0 0 0 1 1 0 1 0;
];

% Number of Monte-Carlo Markov Chains simulation experiments
num_experiments = 5000;

% Number of steps in a single simulation before the experiment ends
max_num_tpe= 200;
min_num_tpe= 10;

% Values for r
R = 0.5:0.5:5;

res = zeros(max_num_tpe - min_num_tpe, 10);

% Run the experiments
i = 1;
for r = R

    j = 1;
    for trials = min_num_tpe:max_num_tpe
        [X, Y] = mcmc(G, r, '01000000', num_experiments, trials);
        res(i, j) = X/num_experiments;
        j = j + 1;
    end
    i = i + 1;

end

hold on
for i = 1:length(R)
    plot(res(i, :))
end
hold off

labels = string(R);
legend(labels)
