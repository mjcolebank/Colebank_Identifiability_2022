% Run a local sensitivity analysis of the triseg model parameterized for
% mouse hemodynamics
clear; clc; close all;
%%
NC     = 30;
data = MI_Sham_data_dynamic;
[IC,pars] = get_pars_newsarc(data);
Tend   = NC.*pars.T;
dt     = 1e-4;

tspace = 0:dt:Tend;

f = @(q) call_model_sens(q,IC,tspace,NC,data);
STEP = 0.01;
NORMALIZE = 0;
par_ids = [1:10 13 16:20 22:23 27 30:34 36:37 38:49];
[sens,nom_sol] = local_sensitivity_log(f,pars,STEP,NORMALIZE,par_ids);

