% Construct upper and lower bounds that are specific to each parameter
% Author: MJ Colebank
% Last edited: 3/10/2022
%
% THIS CODE COMES WITH NO GUARANTEES

function [upper,lower] = get_morris_bounds(pars,ids)
npar = length(pars);
upper = ones(1,npar);
lower = ones(1,npar);

%Pars 1-10 are related to wall volume and reference area
upper(1:5) = 1.5;
lower(1:5) = 0.5;

upper(6:10) = 1.25;
lower(6:10) = 0.75;

% Pars 11-24 are related to atrial sarcomere dynamics
sub_ids = [13 16 17 18 19 20 22 23 24];
upper(sub_ids) = 1.25;
lower(sub_ids) = 0.75;

% Pars 25-37 are related to ventricular sarcomere dynamics
sub_ids = [27 30 31 32 33 34 36 37];
upper(sub_ids) = 1.25;
lower(sub_ids) = 0.75;

% Pars 38-49 are cardiovascular system parameters
% Compartment and valve resistances
sub_ids = [38:45];
upper(sub_ids) = 3.0;
lower(sub_ids) = 0.5;

% Compartment compliances
sub_ids = [46:49];
upper(sub_ids) = 3.0;
lower(sub_ids) = 0.1;

%now return values
upper = upper(ids).*pars(ids);
lower = lower(ids).*pars(ids);
end

%          'Lsref_{v}','Lsiso_{v}','vmax_{v}','Lsc0_{v}','Crest_{v}',...
%           'tauR_{v}','tauD_{v}','tauSC_{v}','sigact_{v}','sigpas_{v}',...
%           'Lsref-pas_{v}','Ls-pas-stiff{v}','k1_{v}', ... 
%           'Ra_val','Rm_val','Rp_val','Rt_val','Rvc','Rpv',...
%           'Rs','Rp','Csa','Csv','Cpa','Cpv'};