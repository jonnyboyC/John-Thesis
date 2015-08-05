% direct{1} = 'D:\thesis\PIVData\mixing layer\Baseline';
direct{1} = 'D:\thesis\PIVData\mixing layer\Forced';

for k = 1
    setup_pod;
    setup_proj;
    
    problem_mixing.direct = direct{k};
    
    
    for j = 1
        problem_mixing.num_clusters = 8 + (j-1)*2;
        POD_Gen(problem_mixing);
    end
    
    problem_mixing.int_time = 7200;
    problem_mixing.tspan = 'test';
    problem_mixing.num_modesG = 4:6;
    Galerkin_Proj(problem_mixing);
    
    problem_mixing.num_modesG = 9:11;
    Galerkin_Proj(problem_mixing);
    
    problem_mixing.num_modesG = 14:16;
    Galerkin_Proj(problem_mixing);
end