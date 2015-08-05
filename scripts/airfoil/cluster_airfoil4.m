addpath(genpath('/home/chabotja/Thesis/John-Thesis'))

direct = '/home/chabotja/Thesis/PIVData/Airfoil/Test_03';

setup_pod;
setup_proj;

problem_airfoil.direct = direct;


% for j = 1
%     problem_airfoil.num_clusters = 8 + (j-1)*2;
%     POD_Gen(problem_airfoil);
% end

problem_airfoil.int_time = 3600;
problem_airfoil.tspan = {'test', 200, 4};
problem_airfoil.num_modesG = 4:6;
Galerkin_Proj(problem_airfoil);
 
problem_airfoil.num_modesG = 9:11;
Galerkin_Proj(problem_airfoil);

problem_airfoil.num_modesG = 14:16;
Galerkin_Proj(problem_airfoil);

problem_airfoil.models = {'GM', 'GM1'};
problem_airfoil.num_modes = 5;
Mod_POD(problem_airfoil);

problem_airfoil.num_modes = 8;
Mod_POD(problem_airfoil);