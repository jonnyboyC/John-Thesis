problem_base1.num_pods = 10;
problem_base1.tspan = 0:0.001:5;
problem_base1.direct = 'D:\thesis\PIVData\Baseline';
problem_base1.Re0_gen = @Re0_gen;

problem_base2.num_pods = 10;
problem_base2.tspan = 0:0.001:5;
problem_base2.direct = 'D:\shear layer\PIVData\Baseline1';
problem_base2.Re0_gen = @Re0_gen;

problem_forced1.num_pods = 10;
problem_forced1.tspan = 0:0.001:5;
problem_forced1.direct = 'D:\thesis\PIVData\Forced';
problem_forced1.Re0_gen = @Re0_gen;

problem_forced2.num_pods = 10;
problem_forced2.tspan = 0:0.001:5;
problem_forced2.direct = 'D:\shear layer\PIVData\Forced1';
problem_forced2.Re0_gen = @Re0_gen;