function likelihood = calc_likelihood(groups, stoch)
% CALC_PROBABILITY calculate the log probability of the markov chain for a
% given transitions matrix
%
% prob = CALC_PROBABILITY(stoch, groups)

transitions = length(groups);

likelihood = 0;
for i = 2:transitions
    likelihood = likelihood+log(stoch(groups(i-1), groups(i)));
end
end