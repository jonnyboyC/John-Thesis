gcp;




%%
load('D:\thesis\PIVData\mixing layer\Forced\POD Data\POD_run_49218.mat');

modal_amp = results_pod.modal_amp(:,2:11);

figure;
eval_mixing_forced1 = evalclusters(modal_amp, 'gmdistribution', 'gap', 'KList', 4:16);
h = plot(eval_mixing_forced1);
title('Mixing Forced Gap');
saveas(h, 'eval_mixing_forced1', 'fig');

figure;
eval_mixing_forced2 = evalclusters(modal_amp, 'gmdistribution', 'silhouette', 'KList', 4:16);
h = plot(eval_mixing_forced2);
title('Mixing Forced silhouette');
saveas(h, 'eval_mixing_forced2', 'fig');

figure;
eval_mixing_forced3 = evalclusters(modal_amp, 'gmdistribution', 'CalinskiHarabasz', 'KList', 4:16);
h = plot(eval_mixing_forced3);
title('Mixing Forced Calinski Harabasz');
saveas(h, 'eval_mixing_forced3', 'fig');

clear results_pod results_clust


%%
load('D:\thesis\PIVData\mixing layer\Baseline\POD Data\POD_run_61768.mat');

modal_amp = results_pod.modal_amp(:,2:11);

figure;
eval_mixing_base1 = evalclusters(modal_amp, 'gmdistribution', 'gap', 'KList', 4:16);
h = plot(eval_mixing_base1);
title('Mixing Base Gap');
saveas(h, 'eval_mixing_base1', 'fig');

figure;
eval_mixing_base2 = evalclusters(modal_amp, 'gmdistribution', 'silhouette', 'KList', 4:16);
h = plot(eval_mixing_base2);
title('Mixing Base silhouette');
saveas(h, 'eval_mixing_base2', 'fig');

figure;
eval_mixing_base3 = evalclusters(modal_amp, 'gmdistribution', 'CalinskiHarabasz', 'KList', 4:16);
h = plot(eval_mixing_base3);
title('Mixing Base Calinski Harabasz');
saveas(h, 'eval_mixing_base3', 'fig');

clear results_pod results_clust

%%

load('D:\thesis\PIVData\airfoil\Test_2\POD Data\POD_run_34400.mat');

modal_amp = results_pod.modal_amp(:,2:11);

figure;
eval_airfoil_base1 = evalclusters(modal_amp, 'gmdistribution', 'gap', 'KList', 4:16);
h = plot(eval_airfoil_base1);
title('Arifoil Base Gap');
saveas(h, 'eval_airfoil_base1', 'fig');

figure;
eval_airfoil_base2 = evalclusters(modal_amp, 'gmdistribution', 'silhouette', 'KList', 4:16);
h = plot(eval_airfoil_base2);
title('Arifoil Base silhouette');
saveas(h, 'eval_airfoil_base2', 'fig');

figure;
eval_airfoil_base3 = evalclusters(modal_amp, 'gmdistribution', 'CalinskiHarabasz', 'KList', 4:16);
h = plot(eval_airfoil_base3);
title('Arifoil Base Calinski Harabasz');
saveas(h, 'eval_airfoil_base3', 'fig');

clear results_pod results_clust


%%
load('D:\thesis\PIVData\airfoil\Test_1\POD Data\POD_run_96670.mat');

modal_amp = results_pod.modal_amp(:,2:11);

figure;
eval_airfoil_forced1 = evalclusters(modal_amp, 'gmdistribution', 'gap', 'KList', 4:16);
h = plot(eval_airfoil_forced1);
title('Arifoil Forced Gap');
saveas(h, 'eval_airfoil_forced1', 'fig');

figure;
eval_airfoil_forced2 = evalclusters(modal_amp, 'gmdistribution', 'silhouette', 'KList', 4:16);
h = plot(eval_airfoil_forced2);
title('Arifoil Forced silhouette');
saveas(h, 'eval_airfoil_forced2', 'fig');

figure;
eval_airfoil_forced3 = evalclusters(modal_amp, 'gmdistribution', 'CalinskiHarabasz', 'KList', 4:16);
h = plot(eval_airfoil_forced3);
title('Arifoil Forced Calinski Harabasz');
saveas(h, 'eval_airfoil_forced3', 'fig');

clear results_pod results_clust

%%
load('D:\thesis\PIVData\cavity\2005_09_16\M030f0000v000a\POD Data\POD_run_49728');

modal_amp = results_pod.modal_amp(:,2:11);

figure;
eval_cavity_base1 = evalclusters(modal_amp, 'gmdistribution', 'gap', 'KList', 4:16);
h = plot(eval_cavity_base1);
title('Cavity Base Gap');
saveas(h, 'eval_cavity_base1', 'fig');

figure;
eval_cavity_base2 = evalclusters(modal_amp, 'gmdistribution', 'silhouette', 'KList', 4:16);
h = plot(eval_cavity_base2);
title('Cavity Base silhouette');
saveas(h, 'eval_cavity_base2', 'fig');

figure;
eval_cavity_base3 = evalclusters(modal_amp, 'gmdistribution', 'CalinskiHarabasz', 'KList', 4:16);
h = plot(eval_cavity_base3);
title('Cavity Base Calinski Harabasz');
saveas(h, 'eval_cavity_base3', 'fig');

clear results_pod results_clust

%%
load('D:\thesis\PIVData\cavity\2005_09_16\M030f1830v400\POD Data\POD_run_36109');

modal_amp = results_pod.modal_amp(:,2:11);

figure;
eval_cavity_forced1 = evalclusters(modal_amp, 'gmdistribution', 'gap', 'KList', 4:16);
h = plot(eval_cavity_forced1);
title('Cavity Forced Gap');
saveas(h, 'eval_cavity_forced1', 'fig');

figure;
eval_cavity_forced2 = evalclusters(modal_amp, 'gmdistribution', 'silhouette', 'KList', 4:16);
h = plot(eval_cavity_forced2);
title('Cavity Forced silhouette');
saveas(h, 'eval_cavity_forced2', 'fig');

figure;
eval_cavity_forced3 = evalclusters(modal_amp, 'gmdistribution', 'CalinskiHarabasz', 'KList', 4:16);
h = plot(eval_cavity_base3);
title('Cavity Forced Calinski Harabasz');
saveas(h, 'eval_cavity_forced3', 'fig');

clear results_pod results_clust


%%
load('D:\thesis\PIVData\DNS jet\streamwise\POD Data\POD_run_82860');

modal_amp = results_pod.modal_amp(:,2:11);

figure;
eval_jet_stream1 = evalclusters(modal_amp, 'gmdistribution', 'gap', 'KList', 4:16);
h = plot(eval_jet_stream1);
title('Jet Base Gap');
saveas(h, 'eval_jet_stream1', 'fig');

figure;
eval_jet_stream2 = evalclusters(modal_amp, 'gmdistribution', 'silhouette', 'KList', 4:16);
h = plot(eval_jet_stream2);
title('Jet Base silhouette');
saveas(h, 'eval_jet_stream2', 'fig');

figure;
eval_jet_stream3 = evalclusters(modal_amp, 'gmdistribution', 'CalinskiHarabasz', 'KList', 4:16);
h = plot(eval_jet_stream3);
title('Jet Base Calinski Harabasz');
saveas(h, 'eval_jet_stream3', 'fig');

clear results_pod results_clust

%%
load('D:\thesis\PIVData\DNS jet\cross_stream2\POD Data\POD_run_4518.mat');

modal_amp = results_pod.modal_amp(:,2:11);

figure;
eval_jet_span1 = evalclusters(modal_amp, 'gmdistribution', 'gap', 'KList', 4:16);
h = plot(eval_jet_span1);
title('Jet spanwise-normal plane Gap');
saveas(h, 'eval_jet_span1', 'fig');

figure;
eval_jet_span2 = evalclusters(modal_amp, 'gmdistribution', 'silhouette', 'KList', 4:16);
h = plot(eval_jet_span2);
title('Jet spanwise-normal plane silhouette');
saveas(h, 'eval_jet_span2', 'fig');

figure;
eval_jet_span3 = evalclusters(modal_amp, 'gmdistribution', 'CalinskiHarabasz', 'KList', 4:16);
h = plot(eval_jet_span3);
title('Jet spanwise-normal plane Calinski Harabasz');
saveas(h, 'eval_jet_span3', 'fig');

clear results_pod results_clust

save('num_clusters');


