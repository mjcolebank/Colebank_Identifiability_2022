%% Run profile likelihood to look at identifiable parameters
% Author: MJ Colebank
% Last edited: 3/10/2022
%
% THIS CODE COMES WITH NO GUARANTEES
clear; clc; close all;
%%
% Want to generate "synthetic data" from regions where our parameters give
% us physiological output
data = MI_Sham_data_dynamic;
[IC,pars] = get_pars_newsarc(data);
par_ids = [1:10 13 16:20 22:23 27 30:34 36:37 38:49]; % Parameters we may change
new_ids = [2 4 5 6 7 9 20 21 22 23 33 34 37];         % Parameters to infer
par_ids = par_ids(new_ids);
NC = 30;
dt     = 1e-4;
T = pars.T;
Tend   = NC.*pars.T;
tspace = 0:dt:Tend;

% Change parameters from structure -> vector
temp = [];
temp(1:10)  = pars.V;
temp(11:24) = pars.SarcA;
temp(25:37) = pars.SarcV;
temp(38:49) = pars.CV;
temp(50)    = pars.T;
pars        = temp;


pars0   = pars;
num_par = length(par_ids);
q0      = 0.*pars(par_ids)+1;   %for scaled parameters
sc      = pars(par_ids);
%% Load in nominal fits/synthetic data
data = load('model_nom_e12.mat');
upper_all = q0.*5;
lower_all = q0.*0.05;

nsamp = 50; %Grid for parameter value
likelihood = zeros(length(q0),nsamp);
par_set    = zeros(length(q0),length(q0),nsamp);

options=optimoptions('lsqnonlin', 'Display','final', 'Algorithm','trust-region-reflective', ...
    'StepTolerance',  1.0000e-8, 'OptimalityTolerance', 1.0000e-8, ...
    'MaxFunctionEvaluations', 8000, 'MaxIterations', 8000);

% DECIDE WHICH RESIDUAL TO USE
rflag = 4;
F  = @(q) call_model_lsqnon(q,0,pars0,par_ids,data,IC,[],tspace,rflag);
%%
clc;

% parpool(12); % Can run in parallel if necessary
% parfor i=1:length(q0)
for i=1:length(q0)
    q_space = linspace(lower(i),upper(i),nsamp);
    q_ids = 1:length(q0);
    q_ids = q_ids(q_ids~=i); % All other indices
    upper_i = upper(q_ids);  % Bounds for ~i parameters
    lower_i = lower(q_ids);  %
    [likelihood_i,par_set_i] = call_loops_lsqnl(rflag,q_space,q_ids,pars0,par_ids, upper_all,lower_all,q0,nsamp,data,IC,sc,tspace,i,options);
    likelihood(i,:) = likelihood_i;
    par_set(i,:,:)  = par_set_i;
end
filename = 'likelihood_lsqnonlin_R4';
save(filename);
