load('/home/chabotja/Thesis/PIVData/Cavity/direct2');

addpath(genpath('/home/chabotja/Thesis/John-Thesis'))

setup_pod;
for i = 1:length(direct)
    for j = 1:3
        problem_cavity.num_clusters = 8 + (j-1)*2;
        problem_cavity.direct = direct{i};
        POD_Gen(problem_cavity);
    end
end