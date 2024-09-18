% This script runs a test for solving the Ax=b for the classic graph and
% getting the fixation probabilities for the nodes u and v given various
% values for r

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

% Initialization for the values of r
num_iters = 50;
step = 0.001;
start = 1;

low = start*step;
up = (num_iters+start-1)*step;
R = low:step:up;
res = zeros((length(G)), length(R), 2);

% Run the test
for threshold = 1:(length(G))
    for i = 1:length(R)
        r = R(i);
        [A, b, idxs_dict] = generateA2(G, r, threshold);
        x = mldivide(A, b);

        u = idxs_dict('10000000');
        v = idxs_dict('01000000');

        res(threshold, i, 1) = x(u);
        res(threshold, i, 2) = x(v);
    end
end
