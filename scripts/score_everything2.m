% direct{1} = 'D:\thesis\PIVData\airfoil\Test';
% direct{2} = 'D:\thesis\PIVData\airfoil\Test_2';
% direct{3} = 'D:\thesis\PIVData\airfoil\Test_12';
% direct{4} = 'D:\thesis\PIVData\airfoil\Test_14';
% direct{5} = 'D:\thesis\PIVData\airfoil\Test_24';
% direct{6} = 'D:\thesis\PIVData\airfoil\Test_26';
% direct{7} = 'D:\thesis\PIVData\airfoil\Test_36';
% direct{8} = 'D:\thesis\PIVData\airfoil\Test_38';
% direct{9} = 'D:\thesis\PIVData\airfoil\Test_49';
% direct{10} = 'D:\thesis\PIVData\airfoil\Test_50';
% direct{11} = 'D:\thesis\PIVData\airfoil\Test_51';
% direct{12} = 'D:\thesis\PIVData\airfoil\Test_52';
direct{1} = 'D:\thesis\PIVData\airfoil\Test_1';
direct{2} = 'D:\thesis\PIVData\airfoil\Test_13';

run_num(1) = 7929;
run_num(2) = 7929;
run_num(3) = 7929;
run_num(4) = 7929;
run_num(5) = 7929;
run_num(6) = 7929;
run_num(7) = 7929;
run_num(8) = 7929;
run_num(9) = 7929;
run_num(10) = 7929;
run_num(11) = 7929;
run_num(12) = 7929;

target_freq = [1250, 1250];
num_clusters = [4, 8, 10, 12];

setup_score;
problem_score.num_cores = 1;
problem_score.target_freq = 0;
problem_score.score_mod = true;

for j = 1:2
    for k = 1:2
        if k == 1
            problem_score.steady = true;
        else
            problem_score.steady = false;
        end
        for i = 1:4
            problem_score.num_clusters =  num_clusters(i);
            problem_score.direct = direct{j};
            problem_score.target_freq = target_freq(j);
            problem_score.run_num = run_num(j);
            score_all(problem_score);
            
        end
    end
end