setup_pod;
problem_airfoil.load_raw = true;
problem_airfoil.direct = 'D:\thesis\PIVData\airfoil\Test_1';
problem_jet.load_raw = true;
problem_jet.direct = 'D:\thesis\PIVData\DNS jet\streamwise';
problem_mixing.load_raw = true;
problem_mixing.direct = 'D:\thesis\PIVData\mixing layer\Baseline';
problem_cavity.load_raw = true;
problem_cavity.direct = 'D:\thesis\PIVData\cavity\2005_09_16\M030f0000v000a';


for i = 1:4
    if i == 1
        delete(gcp);
    else
        parpool('local', i-1);
    end
    
    problem_mixing.num_cores = i;
    tic;
    POD_Gen(problem_mixing);
    time.mixing(i) = toc;
    
    problem_mixing.num_cores = i;
    tic;
    POD_Gen(problem_cavity);
    time.cavity(i) = toc;
    
    problem_mixing.num_cores = i;
    tic;
    POD_Gen(problem_jet);
    time.jet(i) = toc;
    
    problem_mixing.num_cores = i;
    tic;
    POD_Gen(problem_airfoil);
    time.airfoil(i) = toc;
end