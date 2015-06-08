problem_mixing.num_images = 2000;
problem_mixing.load_raw = false;
problem_mixing.save_pod = true;
problem_mixing.l_scale = 0.3048;
problem_mixing.xy_units = 'mm';
problem_mixing.non_dim = false;
problem_mixing.u_scale_gen = @u_scale_gen_shear;
problem_mixing.save_figures = {'fig', 'jpg'};
problem_mixing.flip = [true, false, true, false];

problem_airfoil.num_images = 2000;
problem_airfoil.load_raw = false;
problem_airfoil.save_pod = true;
problem_airfoil.save_figures = {'jpg'};
problem_airfoil.l_scale = 0.2032;
problem_airfoil.u_scale_gen = @u_scale_gen_airfoil;
problem_airfoil.flip = [true, false, true, false];

problem_cavity.num_images = 2000;
problem_cavity.load_raw = false;
problem_cavity.save_pod = true;
problem_cavity.save_figures = {'jpg'};
problem_cavity.u_scale_gen = @u_scale_gen_shear;
problem_cavity.l_scale = 0.0127;

