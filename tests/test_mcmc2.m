% This script runs a test that finds the solutions for the full Ax=b
% problem and runs Monte-Carlo Markov Chains simulations. The purpose of
% this script is to perform a comparison at the results given by both
% methods

% Initial matrix for the graph
G_sample = [
0 1 1 1 0 0 0 0;
1 0 1 0 0 1 0 0;
1 1 0 1 0 0 0 0;
1 0 1 0 0 0 0 1;
0 0 0 0 0 1 1 1;
0 1 0 0 1 0 1 0;
0 0 0 0 1 1 0 1;
0 0 0 1 1 0 1 0;
];

% Final matrix for the problem. It has the initial matrix 4 times
G_final = [G_sample, G_sample; G_sample, G_sample];
n_max = 14;
r = 1.5;

for n = 2:n_max
    disp("----------------------------------------")
    fprintf("G size: %d\n", n);
    G = G_final(1:n, 1:n);

    % Test 1: Simple solution
    fprintf("Simple solution: ");
    [A, b] = generateA(G, r);
    x = mldivide(A, b);
    disp("Done")

    % Test 2: Simple solution with thresholds
    for t = 1:(length(G)-1)
        fprintf("Simple solution with threshold=%d: ", t)
        [A, b] = generateA(G, r, t);
        x = mldivide(A, b);
        disp("Done")
    end

    % Test 3: Monte Carlo Markov Chain
    fprintf("Monte Carlo Markov Chain: ")
    str = pad('01', n, '0');
    [X, Y] = mcmc(G, r, str, 100, 100);
    disp("Done");
end
