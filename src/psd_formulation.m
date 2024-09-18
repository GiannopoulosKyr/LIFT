% This script contains a function to perform Positive Semi-Definite programming
% with various Lassere constraints for getting a solution to the problem 
% Pg=d, where P and d are the values returned by executing generateA.m
% with extinction set to 1.
%
% Arguments for this function are the matrix P and vector d, the
% representations of the various node subsets and the Lassere levels l
%
% Returns the solution to this problem and the matrix X for initializing 
% PSD problem with Lassere.

function [g, X] = psd_formulation2(P, d, reprs, l)
    %%& Initialize
    N = size(P, 1);

    %%% Lassere condition pairs initialization
    % Keep all the combinations that have cardinality of mutants less than
    % or equal to l
    appr = reprs(sum(reprs == '1', 2) <= l, :);

    % Create all the pairs that we will study
    n_units = size(appr, 1);
    pairs1 = nchoosek(1:n_units, 2);    % All pairs, where the elements are not the same
    pairs2 = zeros(n_units, 2);         % All pairs, where the elements are the same 

    for i=1:n_units
        pairs2(i,:) = [i, i];
    end

    %%% Calculate the unions
    % Convert binary representations to decimal
    decimal_appr = bin2dec(appr);

    % Initialize unions vector
    unions1 = zeros(size(pairs1, 1), 1);
    unions2 = zeros(n_units, 1);

    % For each pair, calculate the intersection (the nodes that are mutants
    % on both combinations
    for i=1:size(pairs1, 1)
        unions1(i) = bitor(decimal_appr(pairs1(i, 1)), decimal_appr(pairs1(i, 2)));
    end

    for i=(1:(size(pairs2, 1)))
        unions2(i) = bitor(decimal_appr(pairs2(i, 1)), decimal_appr(pairs2(i, 2)));
    end

    % Find the unique interections and store the idxs of the unions
    % in a vector
    unique_unions = unique([unions1; unions2]);

    %%% Begin Semi-Definite Programming
    cvx_begin sdp
        variable X(N, N)    % Initialize the target matrix
        g = diag(X);        % Store the diagonal of this matrix in a separate variable
        % Initialize a vector that contains all the possible values the
        % elemensts of X can take
        N_vals = size(unique_unions, 1);
        variable vals(N_vals, 1)

        minimize(sum(diag(X)))  % Target

        % Conditions - Part 1
        P * g == d;
        diag(X) >= zeros(N, 1);
        diag(X) <= ones(N, 1);
        X(1,1) == 1;

        % Conditions - Part 2 (Lassere)
        for i = 1:size(unique_unions, 1)
            % For each unique interesection, find the pairs that result in
            % this intersection and they must all have the same value
            % (which means they must have a single value from the values
            % vector, which contains all the possible values X can take)

            % Non-identical pairs
            tmp = pairs1(unions1 == unique_unions(i), :);
            linear_idxs = sub2ind(size(X), tmp(:, 1), tmp(:, 2));
            X(linear_idxs) == vals(i)

            % Non-identical pairs flipped
            tmp = fliplr(tmp);
            linear_idxs = sub2ind(size(X), tmp(:, 1), tmp(:, 2));
            X(linear_idxs) == vals(i)

            % Identical pairs
            tmp = pairs2(unions2 == unique_unions(i), :);
            linear_idxs = sub2ind(size(X), tmp(:, 1), tmp(:, 2));
            X(linear_idxs) == vals(i)
            

            % tmp_linear_idxs = linear_indices(tmp);
            % 
            % tmp_linear_idxs = linear_indices(unions == unique_unions(i));
            % X(tmp_linear_idxs) == vals(i);
        end

        % SDP
        X >= 0;
    cvx_end
end
