direct{1} = 'D:\thesis\PIVData\DNS jet\streamwise';
direct{2} = 'D:\thesis\PIVData\DNS jet\cross_stream1';
direct{3} = 'D:\thesis\PIVData\DNS jet\cross_stream2';
direct{4} = 'D:\thesis\PIVData\DNS jet\cross_stream3';

for k = 1:4
    setup_pod;
    setup_proj;
    
    problem_jet.load_raw = true;
    problem_jet.update_bnds = true;
    problem_jet.cluster = false;
    problem_jet.direct = direct{k};
    
    
    for j = 1
        problem_jet.num_clusters = 8 + (j-1)*2;
        POD_Gen(problem_jet);
    end
    
    problem_jet.classify_sim = false;
    problem_jet.int_time = 3600;
    problem_jet.tspan = {'test'};
    problem_jet.num_modesG = 4:6;
    Galerkin_Proj(problem_jet);
    
    problem_jet.num_modesG = 9:11;
    Galerkin_Proj(problem_jet);
    
    problem_jet.num_modesG = 14:16;
    Galerkin_Proj(problem_jet);
end