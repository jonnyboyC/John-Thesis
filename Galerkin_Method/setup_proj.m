problem_base1.num_modesG = 10;
problem_base1.tspan = 0:0.0001:1;
problem_base1.direct = 'D:\thesis\PIVData\shear layer\Baseline';
problem_base1.Re0_gen = @Re0_gen_shear;

problem_base2.num_modesG = 10;
problem_base2.tspan = 0:0.0001:1;
problem_base2.direct = 'D:\shear layer\PIVData\shear layer\Baseline1';
problem_base2.Re0_gen = @Re0_gen_shear;

problem_forced1.num_modesG = 10;
problem_forced1.tspan = 0:0.0001:1;
problem_forced1.direct = 'D:\thesis\PIVData\shear layer\Forced';
problem_forced1.Re0_gen = @Re0_gen_shear;

problem_forced2.num_modesG = 10;
problem_forced2.tspan = 0:0.0001:1;
problem_forced2.direct = 'D:\shear layer\PIVData\shear layer\Forced1';
problem_forced2.Re0_gen = @Re0_gen_shear;

problem_airfoil.num_modesG = 10;
problem_airfoil.tspan = 0:0.0001:1;
problem_airfoil.Re0_gen = @Re0_gen_airfoil;
