%% Morris screening
% Author: MJ Colebank
% Last edited: 3/10/2022
%
% THIS CODE COMES WITH NO GUARANTEES

clc; clear; 
%%
Names = {'VW_{la}','VW_{lv}','VW_{ra}','VW_{rv}','VW_{s}',...
          'Amref_{la}','Amref_{lv}','Amref_{ra}','Amref_{rv}','Amref_{s}',...
          'Lsref_{a}','Lsiso_{a}','vmax_{a}','Lsc0_{a}','Crest_{a}',...
          'tauR_{a}','tauD_{a}','tauSC_{a}','sigact_{a}','sigpas_{a}',...
          'Lsref-pas_{a}','Ls-pas-stiff{a}','k1_{a}','toffset',...
         'Lsref_{v}','Lsiso_{v}','vmax_{v}','Lsc0_{v}','Crest_{v}',...
          'tauR_{v}','tauD_{v}','tauSC_{v}','sigact_{v}','sigpas_{v}',...
          'Lsref-pas_{v}','Ls-pas-stiff{v}','k1_{v}', ... 
          'Ra_val','Rm_val','Rp_val','Rt_val','Rvc','Rpv',...
          'Rs','Rp','Csa','Csv','Cpa','Cpv'};
%%

%% Initialize things for the model
NC     = 60;
data = MI_Sham_data_dynamic;
[IC,pars] = get_pars_newsarc(data);
Tend   = NC.*pars.T;

dt     = 1e-4;
T = pars.T;
tspace = 0:dt:Tend;
% Change parameters from structure -> vector
temp = [];
temp(1:10) = pars.V;
temp(11:24) = pars.SarcA;
temp(25:37) = pars.SarcV;
temp(38:49) = pars.CV;
temp(50)    = pars.T;
pars = temp;
par_ids = [1:10 13 16:20 22:23 27 30:34 36:37 38:49];
par_nom = pars;
f = @(q) call_model_Morris(q,IC,tspace,par_ids,par_nom);
%%
q_nom = pars(par_ids);
r = 120;               % Max samples until we reach r
smallR = 100;           % Number of samples we want
p = length(par_ids);   % Number of parameters
l = 60;                % Number of levels
delta = l./(2*(l-1));  % Step Size
% MJC 9/13/2021 - new file to get more parameter specific bounds:
[upper,lower] = get_morris_bounds(par_nom,par_ids);

d = cell(r,p);
%% Try to use the randomization algorithm
% Note that all parameters are scaled to be in the range 0,1 and then
% rescaled in the model evaluation.
A = zeros(p+1,p);
for i=1:p
    A(i+1:p+1,i) = ones((p-(i-1)),1);
end
X = zeros((p+1).*r,p.*r);
%% 
F_storage = cell(p+1,smallR);
qstar = unifrnd(0,1,r,p);
Jp = ones(p+1,p);
J1 = ones(p+1,1);
P = eye(p,p);
UL_MAT = eye(p).*(upper-lower);
i=1;
func_evals = 1;
while i<=smallR
    exitflag = 0;
    if func_evals == r
        error('Parameter Space Infeasible: Exiting.\n');
    end
    qcurr = qstar(i,:);
    pm1 = rand(p,1);
    Dstar = eye(p).*(pm1 > 0.5) - eye(p).*(pm1 <= 0.5);
    [who,where] = sort(rand(p,1));
    Pstar = P(where,:);
    Astar = J1*qcurr + (delta./2).*(( (2.*A - Jp)*Dstar + Jp))*Pstar;
    C = J1*(lower) + Astar*UL_MAT;
    fpast = f(C(1,:));
    F_storage{1,i} = fpast;
    for j=1:p
        disp([j+1 Names(par_ids(where(j)))])
        fpresent = f(C(j+1,:));
        if isempty(fpresent) || isempty(fpast)
           for s=j:-1:1
              d{i,where(s)} = {}; %Clear all previous entries 
           end
           exitflag = 1;
           break;
        end
        d{i,where(j)} = (fpresent - fpast)./delta;

        fpast = fpresent;
        F_storage{where(j)+1,i} = fpast;
    end
    if exitflag == 1
        exitflag = 0;
    else
        i = i+1;
    end
end
save('Morris_results');

