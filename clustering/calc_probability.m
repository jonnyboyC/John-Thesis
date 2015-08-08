function prob = calc_probability(stoch, groups)
% CALC_PROBABILITY calculate the log probability of the markov chain for a
% given transitions matrix
%
% prob = CALC_PROBABILITY(stoch, groups)
transitions = length(groups);

prob = 0;
for i = 2:transitions
    prob = prob+log(stoch(groups(i-1), groups(i)));
end
end