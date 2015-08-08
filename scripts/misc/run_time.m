setup_pod;
problem_airfoil.load_raw = true;
problem_airfoil.direct = 'D:\thesis\PIVData\airfoil\Test_1';
problem_jet.load_raw = true;
problem_jet.direct = 'D:\thesis\PIVData\DNS jet\streamwise';
problem_mixing.load_raw = true;
problem_mixing.direct = 'D:\thesis\PIVData\mixing layer\Baseline';
problem_cavity.load_raw = true;
problem_cavity.direct = 'D:\thesis\PIVData\cavity\2005_09_16\M030f0000v000a';

% 
% for i = 1:4
%     
% %     problem_mixing.num_cores = i;
% %     tstart = tic;
% %     POD_Gen(problem_mixing);
% %     time_pod.mixing(i) = toc(tstart);
%     
%     problem_cavity.num_cores = i;
%     tic;
%     POD_Gen(problem_cavity);
%     time_pod.cavity(i) = toc;
%     
% %     problem_jet.num_cores = i;
% %     tstart = tic;
% %     POD_Gen(problem_jet);
% %     time_pod.jet(i) = toc(tstart);
%     
% %     problem_airfoil.num_cores = i;
% %     tstart = tic;
% %     POD_Gen(problem_airfoil);
% %     time_pod.airfoil(i) = toc(tstart);
% end

% for i = 1:4
%     
%     problem_mixing.num_cores = i;
%     problem_mixing.num_modesG = 10;
%     problem_mixing.tspan = {'test', 20};
%     problem_mixing.override_coef = true;
%     tstart = tic;
%     Galerkin_Proj(problem_mixing);
%     time_gal.mixing(i) = toc(tstart);
%     
%     problem_cavity.num_cores = i;
%     problem_cavity.num_modesG = 10;
%     problem_cavity.tspan = {'test', 400};
%     problem_cavity.override_coef = true;
%     tic;
%     Galerkin_Proj(problem_cavity);
%     time_gal.cavity(i) = toc;
%     
%     problem_jet.num_cores = i;
%     problem_jet.num_modesG = 10;
%     problem_jet.tspan = {'test'};
%     problem_jet.override_coef = true;
%     tstart = tic;
%     Galerkin_Proj(problem_jet);
%     time_gal.jet(i) = toc(tstart);
%     
%     problem_airfoil.num_cores = i;
%     problem_airfoil.num_modesG = 10;
%     problem_airfoil.tspan = {'test', 20};
%     problem_airfoil.override_coef = true;
% 
%     tstart = tic;
%     Galerkin_Proj(problem_airfoil);
%     time_gal.airfoil(i) = toc(tstart);
% end

setup_mod
for i = 1:4
    
%     problem_mixing.num_cores = i;
%     problem_mixing.num_modes = 5;
%     tstart = tic;
%     Mod_POD(problem_mixing);
%     time_mod.mixing(i) = toc(tstart);
    
    problem_cavity.num_cores = i;
    problem_cavity.num_modes = 5;
    tic;
    Mod_POD(problem_cavity);
    time_mod.cavity(i) = toc;
    
%     problem_jet.num_cores = i;
%     problem_jet.num_modes = 5;
%     tstart = tic;
%     Mod_POD(problem_jet);
%     time_mod.jet(i) = toc(tstart);
    
%     problem_airfoil.num_cores = i;
%     problem_airfoil.num_modes = 5;
%     tstart = tic;
%     Mod_POD(problem_airfoil);
%     time_mod.airfoil(i) = toc(tstart);
end