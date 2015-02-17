function list = list_check(list, correct, var)
% Small function to condense the setup defaults files
correct_members = ismember(list, correct);
for i = 1:size(correct_members,2)
    if ~correct_members(i)
        fprintf('%s is not a correct input for problem.%d\n', list{i}, var);
    end
end
list = list(correct_members);
end