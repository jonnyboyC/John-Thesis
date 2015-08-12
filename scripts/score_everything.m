

problem_score.num_cores = 1;
for i = 1:3
    problem_score.num_clusters = 6 + i*2;
    
    problem_score.direct = 'D:\thesis\PIVData\mixing layer\Forced';
    problem_score.target_freq = 60;
    score_all(problem_score);
    
    problem_score.direct = 'D:\thesis\PIVData\mixing layer\Baseline';
    problem_score.target_freq = 0;
    score_all(problem_score);
end


direct{1} = 'D:\thesis\PIVData\DNS jet\streamwise';
direct{2} = 'D:\thesis\PIVData\DNS jet\cross_stream1';
direct{3} = 'D:\thesis\PIVData\DNS jet\cross_stream2';
direct{4} = 'D:\thesis\PIVData\DNS jet\cross_stream3';

for i = 1:3
    for j = 1:4
        problem_score.direct = direct{j};
        problem_score_target_freq = 0;
        score_all(problem_score);
    end
end




