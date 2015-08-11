direct{1} = 'D:\thesis\PIVData\DNS jet\streamwise';
direct{2} = 'D:\thesis\PIVData\DNS jet\cross_stream1';
direct{3} = 'D:\thesis\PIVData\DNS jet\cross_stream2';
direct{4} = 'D:\thesis\PIVData\DNS jet\cross_stream3';

for k = 1:4
    setup_pod;
    setup_proj;
    
    problem_jet.update_bnds = false;
    problem_jet.load_raw = false;
    problem_jet.direct = direct{k};
    

    problem_jet.num_clusters = 10;
    POD_Gen(problem_jet);

   
    problem_jet.int_time = 3600;
    problem_jet.tspan = {'test', 1, 4};
    problem_jet.num_modesG = 4:6;
    Galerkin_Proj(problem_jet);
    
    problem_jet.num_modesG = 9:11;
    Galerkin_Proj(problem_jet);
    
    problem_jet.num_modesG = 14:16;
    Galerkin_Proj(problem_jet);
    
    problem_jet.models = {'GM', 'GM1'};
    problem_jet.num_modes = 5;
    Mod_POD(problem_jet);

    problem_jet.num_modes = 8;
    Mod_POD(problem_jet);
end