addpath(genpath('/home/chabotja/Thesis/John-Thesis'))

direct = '/home/chabotja/Thesis/PIVData/Cavity/2005_09_16/M030f0000v000b';

setup_pod;
setup_proj;

problem_cavity.direct = direct;


for j = 1
    problem_cavity.num_clusters = 8 + (j-1)*2;
    POD_Gen(problem_cavity);
end

problem_cavity.int_time = 14400;
problem_cavity.tspan = 'test';
problem_cavity.num_modesG = 4:6;
Galerkin_Proj(problem_cavity);

problem_cavity.num_modesG = 9:11;
Galerkin_Proj(problem_cavity);

problem_cavity.num_modesG = 14:16;
Galerkin_Proj(problem_cavity);