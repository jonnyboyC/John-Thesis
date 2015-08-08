direct{1} = 'D:\thesis\PIVData\mixing layer\Baseline';
direct{2} = 'D:\thesis\PIVData\mixing layer\Forced';

for k = 1:2
    setup_pod;
    setup_proj;
    
    problem_mixing.direct = direct{k};
    
    for j = 1
        problem_mixing.num_clusters = 8 + (j-1)*2;
        POD_Gen(problem_mixing);
    end
    
    problem_mixing.int_time = 3600;
    problem_mixing.tspan = {'test', 400, 4};
    problem_mixing.num_modesG = 4:6;
    Galerkin_Proj(problem_mixing);
    
    problem_mixing.num_modesG = 9:11;
    Galerkin_Proj(problem_mixing);
    
    problem_mixing.num_modesG = 14:16;
    Galerkin_Proj(problem_mixing);
    
    problem_mixing.models = {'GM', 'GM1'};
    problem_mixing.num_modes = 5;
    Mod_POD(problem_mixing);

    problem_mixing.num_modes = 8;
    Mod_POD(problem_mixing);
end