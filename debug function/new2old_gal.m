%% This script simply changes the names of some variables in the workspace 
% so they comply with the naming convenction of previous code development.
clc
clear all

load('D:\shear layer\PIVData\Forced1\POD Data\POD.mat');
muz = {mean_u};
mvz = {mean_v};
phiutc = {pod_u1};
phivtc = {pod_v1};
vt = {vol_frac};
bi = bnd_idx;
sgtc2 = lambda2;
psitc = eig_func;
x = {x};
y = {y};
sa = {dimensions};
t0 = [1 2];
Nsamples = 2000;
N = 42504;

clear mean_u mean_v pod_u1 pod_v1 vol_frac bnd_idx lambda2 eig_func dimensions
