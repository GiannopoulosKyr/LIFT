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

assert(issymmetric(G));
assert(all(diag(G) == 0));

num = length(G);
b = zeros(2^num, 1);
b(2^num) = 1;

num_iters = 10000;
step = 0.0005;
start = 1;

low = start*step;
up = (num_iters+start-1)*step;

R = low:step:up;

res = zeros(length(R), 2);

for i = 1:length(R)

    r = R(i);
    A = generateA(G, r);
    x = A\b;    

    u = bin2dec('0001')+1;
    v = bin2dec('0010')+1;

    res(i, 1) = x(u);
    res(i, 2) = x(v);

end

% Plot results

tiledlayout(3,1)

% Plot1
nexttile
plot(R, res)
legend("u", "v")
title("Plot 1: Fixation probabilities")
xlabel('r')
ylabel('Fixation probabilities')

% Plot2
nexttile
frac = res(:,1) ./ res(:, 2);
plot(R, frac)
y1 = yline(1, '-.');
legend("Ratio")
title("Plot 2: Ratio of fixation probabilities")
xlabel('r')
ylabel('Ratio')

% Plot 3
nexttile
d = res(:,1) - res(:,2);
[max_diff, max_ind] = max(d);
plot(R, d)
y2 = yline(0, '-.');
legend("Difference")
title("Difference of fixation probabilities")
xlabel('Relative fitness')
ylabel('')

max_r = R(max_ind);

fprintf("Maximum r: %f\n", max_r);
