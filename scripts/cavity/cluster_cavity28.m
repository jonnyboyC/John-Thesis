addpath(genpath('/home/chabotja/Thesis/John-Thesis'))

direct = '/home/chabotja/Thesis/PIVData/Cavity/2007_5_4/M030fM0Bcomp';

setup_pod;
setup_proj;

problem_cavity.direct = direct;

problem_cavity.num_clusters = 10;
POD_Gen(problem_cavity);

problem_cavity.int_time = 3600;
problem_cavity.tspan = {'test', 200, 4};
problem_cavity.num_modesG = 4:6;
Galerkin_Proj(problem_cavity);

problem_cavity.num_modesG = 9:11;
Galerkin_Proj(problem_cavity);

problem_cavity.num_modesG = 14:16;
Galerkin_Proj(problem_cavity);

problem_cavity.models = {'GM', 'GM1'};
problem_cavity.num_modes = 5;
Mod_POD(problem_cavity);

problem_cavity.num_modes = 8;
Mod_POD(problem_cavity);