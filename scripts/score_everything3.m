direct{1} = 'D:\thesis\PIVData\cavity\2005_09_14\M030f0000v000a';
direct{2} = 'D:\thesis\PIVData\cavity\2005_09_15\M030f0000v000a';
direct{3} = 'D:\thesis\PIVData\cavity\2005_09_15\M030f0000v000b';
direct{4} = 'D:\thesis\PIVData\cavity\2005_09_15\M030f0000v000c';
direct{5} = 'D:\thesis\PIVData\cavity\2005_09_16\M030f0000v000a';
direct{6} = 'D:\thesis\PIVData\cavity\2005_09_16\M030f0000v000b';
direct{7} = 'D:\thesis\PIVData\cavity\2005_09_19\M030f0000v000';
direct{8} = 'D:\thesis\PIVData\cavity\2005_09_20\M030f1830v400a';
direct{9} = 'D:\thesis\PIVData\cavity\2005_09_20\M030f1830v400b';
direct{10} = 'D:\thesis\PIVData\cavity\2005_09_20\M030f3000v300';
direct{11} = 'D:\thesis\PIVData\cavity\2005_10_14\M030f3000v300';
direct{12} = 'D:\thesis\PIVData\cavity\2005_10_18\M030f1610v250';
direct{13} = 'D:\thesis\PIVData\cavity\2005_10_18\M030f3920v250';
direct{14} = 'D:\thesis\PIVData\cavity\2005_10_20\M030f0000v000a';
direct{15} = 'D:\thesis\PIVData\cavity\2005_10_20\M030f0000v000b';
direct{16} = 'D:\thesis\PIVData\cavity\2005_10_20\M030f2900v600b';
direct{17} = 'D:\thesis\PIVData\cavity\2005_10_20\M030f3000v200';
direct{18} = 'D:\thesis\PIVData\cavity\2005_10_21\M030f1610v170';
direct{19} = 'D:\thesis\PIVData\cavity\2005_10_21\M030f2800v750';
direct{20} = 'D:\thesis\PIVData\cavity\2005_10_21\M030f3000v150';
direct{21} = 'D:\thesis\PIVData\cavity\2005_10_21\M030f3000v200';
direct{22} = 'D:\thesis\PIVData\cavity\2007_5_4\M030fM1BF4';
direct{23} = 'D:\thesis\PIVData\cavity\2007_5_4\M030fM1BF4comp';
direct{24} = 'D:\thesis\PIVData\cavity\2007_5_4\M030fWNcomp3';

run_num{1} = 62695;
run_num{2} = 'first';
run_num{3} = 'first';
run_num{4} = 'first';
run_num{5} = 6400;
run_num{6} = 'first';
run_num{7} = 'first';
run_num{8} = 'first';
run_num{9} = 6400;
run_num{10} = 'first';
run_num{11} = 'first';
run_Num{12} = 'first';
run_num{13} = 'first';
run_num{14} = 'first';
run_num{15} = 'first';
run_num{16} = 'first';
run_num{17} = 'first';
run_num{18} = 'first';
run_num{19} = 'first';
run_num{20} = 'first';
run_num{21} = 'first';
run_num{22} = 'first';
run_num{23} = 'first';
run_num{24} = 'first';

% run_num(10) = 7929;
% run_num(11) = 7929;
% run_num(12) = 7929;

setup_score;
problem_score.num_cores = 2;
problem_score.target_freq = 0;
problem_score.score_mod = false;

for i = 1:3
    problem_score.num_clusters = 6 + 2*i;
    for j = 1:11
        problem_score.direct = direct{j};
        problem_score.run_num = run_num{j};
        score_all(problem_score);
    end
end