problem_mixing.num_images = 2000;
problem_mixing.load_raw = false;
problem_mixing.save_pod = true;
problem_mixing.l_scale = 0.3048;
problem_mixing.xy_units = 'mm';
problem_mixing.non_dim = false;
problem_mixing.u_scale_gen = @u_scale_gen_shear;
problem_mixing.save_figures = {'fig', 'jpg'};
problem_mixing.load_handle = @load_LaVision;
problem_mixing.flip = [true, false, true, false];

problem_airfoil.num_images = 2000;
problem_airfoil.load_raw = false;
problem_airfoil.save_pod = true;
problem_airfoil.save_figures = {'fig', 'jpg'};
problem_airfoil.load_handle = @load_LaVision;
problem_airfoil.xy_units = 'mm';
problem_airfoil.non_dim = false;
problem_airfoil.l_scale = 0.2032;
problem_airfoil.u_scale_gen = @u_scale_gen_airfoil;
problem_airfoil.flip = [true, false, true, false];

problem_cavity.num_images = 2000;
problem_cavity.load_raw = false;
problem_cavity.load_handle = @load_cavity_dat;
problem_cavity.save_pod = true;
problem_cavity.xy_units = 'mm';
problem_cavity.non_dim = false;
problem_cavity.save_figures = {'fig', 'jpg'};
problem_cavity.u_scale_gen = @u_scale_gen_shear;
problem_cavity.l_scale = 0.0127;

problem_jet.num_images = 860;
problem_jet.load_raw = false;
problem_jet.load_handle = @load_jet_dat;
problem_jet.save_pod = true;
problem_jet.xy_units = 'm';
problem_jet.non_dim = true;
problem_jet.save_figures = {'fig', 'jpg'};
problem_jet.u_scale_gen = 3.280839;
problem_jet.l_scale = 1/0.3048;


