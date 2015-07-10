addpath(genpath('/home/chabotja/Thesis/John-Thesis'))

direct = '/home/chabotja/Thesis/PIVData/Cavity/2007_5_4/M030fWNcomp1';

setup_pod;
setup_proj;

problem_cavity.direct = direct;


for j = 1:3
    problem_cavity.num_clusters = 8 + (j-1)*2;
    POD_Gen(problem_cavity);
end

problem_cavity.tspan = 'test';
problem_cavity.num_modes = 4:6;
Galerkin_Proj(problem_cavity);

problem_cavity.num_modes = 9:11;
Galerkin_Proj(problem_cavity);

problem_cavity.num_modes = 14:16;
Galerkin_Proj(problem_cavity);