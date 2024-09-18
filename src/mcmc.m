function [fixation_counter, non_fixation_counter] = mcmc(G, r, S, num_experiments, num_trials)

    assert(length(G) == length(S))
    assert(issymmetric(G))
    assert(sum(S == '1') ~= 0)

    fixation_counter = 0;
    non_fixation_counter = 0;
    num_nodes = length(G);

    for e = 1:num_experiments
        S_copy = S;

        for t = 1:num_trials
            num_mutants = sum(S_copy == '1');
            num_residents = num_nodes - num_mutants;

            probs = ones(num_nodes, 1);
            probs(S_copy == '1') = r;
            probs = probs / (num_residents + r*num_mutants);

            % Random selection of node
            cum_probs = cumsum(probs);
            rand_num = rand;
            node = find(rand_num < cum_probs, 1);

            if S_copy(node) == '1'
                adjacent_cum = (1:(num_nodes-1))/(num_nodes-1);
                adj_rand = rand;
                adj_node = find(adj_rand < adjacent_cum, 1);
                if adj_node >= node
                    adj_node = adj_node + 1;
                end
                
                S_copy(adj_node) = '1';
            else
                ordered_nodes = (1:num_nodes);
                adjacent_nodes = ordered_nodes(G(node) == 1);

                adjacent_cum = (1:length(adjacent_nodes))/(length(adjacent_nodes));
                adj_rand = rand;
                adj_node = adjacent_nodes(find(adj_rand < adjacent_cum, 1));

                S_copy(adj_node) = '0';

            end

            if sum(S_copy == '1') == num_nodes || sum(S_copy == '0') == num_nodes
                break
            end
        end

        if sum(S_copy == '1') == num_nodes
            fixation_counter = fixation_counter + 1;
        else
            non_fixation_counter = non_fixation_counter + 1;
        end

    end

end
