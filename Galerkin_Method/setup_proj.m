problem_mixing.num_modesG = 10;
problem_mixing.tspan = 0:0.0001:1;
problem_mixing.Re0_gen = @Re0_gen_shear;

problem_cavity.num_modesG = 10;
problem_cavity.tspan = 0:0.0001:1;
problem_cavity.Re0_gen = @Re0_gen_cavity;

problem_airfoil.num_modesG = 10;
problem_airfoil.tspan = 0:0.0001:1;
problem_airfoil.Re0_gen = @Re0_gen_airfoil;

% TODO update to actual value
problem_jet.num_modesG = 10;
problem_jet.tspan = 0:0.0001:1;
problem_jet.Re0_gen = @Re0_gen_shear;
