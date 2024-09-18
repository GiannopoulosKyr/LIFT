% Script that runs PSD programming for the problem with Lassere Constraints

% Initial graph
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

% Values for r
R = [0.375, 0.975, 1.275, 1.575, 1.875];

% The set of graphs to be tested
graphs = {G};

% Run the test
res = test(graphs, R);

% The function for test
function [res] = test(graphs, R)
    res = {};

    for i=1:length(graphs)
        G = graphs{i};
        n = size(G, 1);
        res{i} = zeros([length(R), size(G, 1)-1, size(G, 1)-1, 9]);

        for j=1:length(R)
            r = R(j);

            [A, d, reprs] = generateA(G, r, 0, 1);
            g = mldivide(A, d);

            % Truncation level k
            for k=1:(n-1)
                [P1, d1, inv_reprs1, reprs1] = generateA(G, r, k, 1);
                [P2, d2, inv_reprs2, reprs2] = generateA(G, r, k+1, 1);

                g1 = mldivide(P1, d1);
                S = size(P1,1);
                P1 = P2(1:size(P2, 1), :);

                % Lassere level l
                for l=1:k

                    [g2, X] = psd_formulation(P1, d2, reprs2, l);

                    tmp1 = g(1:S) - g2(1:S);
                    tmp2 = g2(1:S) - g1;
                    tmp3 = g(1:S) - g1;

                    res{i}(j, k, l, 1) = sum(tmp1);
                    res{i}(j, k, l, 2) = sum(tmp2);
                    res{i}(j, k, l, 3) = sum(tmp3);
                    res{i}(j, k, l, 4) = mean(tmp1);
                    res{i}(j, k, l, 5) = mean(tmp2);
                    res{i}(j, k, l, 6) = mean(tmp3);
                    res{i}(j, k, l, 7) = all(tmp1 >= 0);
                    res{i}(j, k, l, 8) = all(tmp2 >= 0);
                    res{i}(j, k, l, 9) = all(tmp3 >= 0);
                    res{i}(j, k, l, 10) = res{i}(j, k, l, 4) >= 0;
                    res{i}(j, k, l, 11) = res{i}(j, k, l, 5) >= 0;
                    res{i}(j, k, l, 12) = res{i}(j, k, l, 6) >= 0;

                    full(X)

                end
            end
        end
    end
end
