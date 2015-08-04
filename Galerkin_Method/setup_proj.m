problem_mixing.num_modesG = 10;
problem_mixing.tspan = {'test', 100};
problem_mixing.odesolver = @ode113;
problem_mixing.Re0_gen = @Re0_gen_shear;

problem_cavity.num_modesG = 10;
problem_cavity.tspan = {'test', 100};
problem_cavity.odesolver = @ode113;
problem_cavity.Re0_gen = @Re0_gen_cavity;

problem_airfoil.num_modesG = 10;
problem_airfoil.tspan = {'test', 200};
problem_airfoil.odesolver = @ode113;
problem_airfoil.Re0_gen = @Re0_gen_airfoil;

% TODO update to actual value
problem_jet.num_modesG = 10;
problem_jet.tspan = {'test'};
problem_jet.odesolver = @ode113;
problem_jet.Re0_gen = @Re0_gen_jet;
