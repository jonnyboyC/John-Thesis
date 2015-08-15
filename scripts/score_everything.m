
% num_clusters = [5, 8, 10, 12];
% 
% problem_score.num_cores = 1;
% for i = 3:4
%     for j = 1:2
%         if j == 1
%             problem_score.steady = true;
%         else
%             problem_score.steady = false;
%         end
%         problem_score.num_clusters = num_clusters(i);
%         
%         problem_score.direct = 'D:\thesis\PIVData\mixing layer\Forced';
%         problem_score.target_freq = 60;
%         problem_score.phase = [1,2];
%         score_all(problem_score);
%         
%         problem_score = rmfield(problem_score, 'phase');
%         problem_score.direct = 'D:\thesis\PIVData\mixing layer\Baseline';
%         problem_score.target_freq = 0;
%         score_all(problem_score);
%     end
% end

% problem_score = rmfield(problem_score, 'target_freq');

direct{1} = 'D:\thesis\PIVData\DNS jet\streamwise';
direct{2} = 'D:\thesis\PIVData\DNS jet\cross_stream1';
direct{3} = 'D:\thesis\PIVData\DNS jet\cross_stream2';
direct{4} = 'D:\thesis\PIVData\DNS jet\cross_stream3';

num_clusters = [5, 8, 10, 16];

for i = 4
    for j = 1:4
        for k = 1:2
            if k == 1
                problem_score.steady = true;
            else
                problem_score.steady = false;
            end
            problem_score.direct = direct{j};
            problem_score.num_clusters = num_clusters(i);
            score_all(problem_score);
        end
    end
end




