% Run the DRAM algorithm on the TriSeg model IN PARALLEL
% Author: MJ Colebank
% Last edited: 3/10/2022
%
% THIS CODE COMES WITH NO GUARANTEES
clear; clc;
addpath('DRAM/'); % NOTE must download the MCMSTAT DRAM folder available online
%%
NC   = 30;
data = MI_Sham_data_dynamic;

[IC,pars] = get_pars_newsarc(data);
Tend   = NC.*pars.T;
dt     = 1e-4;
tspace = 0:dt:Tend;

par_ids = [1:10 13 16:20 22:23 27 30:34 36:37 38:49]; %Parameters we may change
test_ids = [2 4 5 6 7 9 20 21 22 23 33 34 37]; %Morris+SVD
IDS = par_ids(test_ids);
IDS = sort(IDS);

% Load model predictions or noisey data
% data = load('model_nom_e12.mat');
data = load('data_noise_std_P1_V1.mat'); %Noise
[~,qall,UB,LB] = get_bounds(pars,IDS);
rflag = 1;
F = @(q,dummy) opt_wrap(q,IDS,qall,IC,tspace,NC,data,rflag);
Jopt = @(q) F(q,0);
%% Run in parallel
nsamp = 1;
chain_storage   = cell(nsamp,1);
results_storage = cell(nsamp,1);
s2_storage      = cell(nsamp,1);

% parpool(nsamp); % Can run in parallel or not
% parfor i = 1:nsamp
for i=1:nsamp
[q0,~,~,~] = get_bounds(pars,IDS);
[results,chain,s2chain] = call_dram(q0,F,data,IDS,UB,LB,i,rflag);
chain_storage{i}  = chain;
results_storage{i}= results;
s2_storage{i}     = s2chain;
end
save('DRAM_Files');
