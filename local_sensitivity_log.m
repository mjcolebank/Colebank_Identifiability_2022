%% Local_sensitivity.m
% Written by: MJ Colebank
% August 26, 2019
% Updated: August 13, 2021
%
% Computes a local sensitivity approximation using centered finite
% differences. Contains the subfunction call_centered_difference.m
%
% INPUTS:
% f - the function/system of interest; could be solution to a ODE/PDE/DAE.
% This is the quantity that will be used in computing the sensitivity, i.e.
% df/dtheta. This can be changed to be the max,min, etc. but must be
% changed before passing f into the function.
%
% pars - the parameters of the system you wish to perturb.
%
% data - the independent variable (either space or time)
%
% STEP - The stepsize applied for the perturbations in the finite
% difference scheme. Note that this "parameter" should be tuned to give
% enough resolution in computing the finite difference, but not so fine
% that your sensitivity is zero.
%
% NORMALIZE - Computes the sensitivities by normalizing them with respect
% to the nominal output. Needs to be 0 or 1 (1 for normalization).
%
% This code comes with no guarantee.
%
%
%

function [sens,nom_sol] = local_sensitivity_log(f,pars,STEP,NORMALIZE,par_ids)

%Depending on form of y, the sensitivity matrix could be dymodel/dpars or
%dr/dpars where r = ymodel - ydata
if nargin==4
    par_ids = 1:length(pars);
end  
pars = struct2cell(pars);
pars = cell2mat(pars)';
nom_sol = f(pars);
% [sens] = call_forward_difference(pars,f,nom_sol,STEP,NORMALIZE,par_ids);
[sens] = call_centered_difference(pars,f,nom_sol,STEP,NORMALIZE,par_ids);
end


function sens = call_forward_difference(pars,f,nom_sol,STEP,NORMALIZE,par_ids)
P = par_ids;
Pn = length(P);
outs = size(nom_sol,1);   % Number of QoI
N    = size(nom_sol,2);   % Number of time points
par0 = log(pars);
if outs==1
    sens = zeros(N,Pn);
    for i=1:Pn %Note: we only run through the parameters in par_ids
        par_plus  = par0; %
        par_plus(P(i)) = par0(P(i))+STEP; %params + increment,Expoentiate
        f_plus = feval(f,exp(par_plus)); %Solution of ODE at the increment
        sens(:,i) = (f_plus - nom_sol)./STEP; %approximate the derivative
    end
else % for more than one QoI
    sens = zeros(outs,N,Pn);
    for i=1:Pn %parfor
%         disp(i);

        par_plus  = par0;
        par_plus(P(i)) = par0(P(i))+STEP;
        disp([P(i) par0(P(i)) par_plus(P(i))]);
        f_plus = feval(f,exp(par_plus));
            
        % Might need to reorient if necessary 
        for j=1:outs
            sens(j,:,i) = (f_plus(j,:) - nom_sol(j,:))./STEP;
        end
    end
end

% If set to 1, normalize the sensitivities with respect to the magnitude of
% the output (e.g., divide by the solution from the non-perturbed
% parameters
if NORMALIZE == 1
    if outs>1
        for k=1:outs
            sens(k,:,:) = sens(k,:,:)./nom_sol(k,:);
        end
    else
        sens = sens./nom_sol';
    end
end
end

function sens = call_centered_difference(pars,f,nom_sol,STEP,NORMALIZE,par_ids)
P = par_ids;
Pn = length(P);
outs = size(nom_sol,1);   % Number of QoI
N    = size(nom_sol,2);   % Number of time points
par0 = log(pars);
if outs==1
    sens = zeros(N,Pn);
    for i=1:Pn %Note: we only run through the parameters in par_ids
        par_plus  = par0;
        par_plus(P(i)) = par0(P(i))+STEP; %params + increment
        par_minus  = par0;
        par_minus(P(i)) = par0(P(i))-STEP; %params + increment
        f_plus = feval(f,exp(par_plus)); %Solution of ODE at the increment
        f_minus = feval(f,exp(par_minus)); %Solution of ODE at the increment
        sens(:,i) = (f_plus - f_minus)./(2.*STEP); %approximate the derivative
    end
else % for more than one QoI
    sens = zeros(outs,N,Pn);
    for i=1:Pn %parfor
        disp(i);

        par_plus  = par0;
        par_plus(P(i)) = par0(P(i))+STEP; %params + increment
        par_minus  = par0;
        par_minus(P(i)) = par0(P(i))-STEP; %params + increment
        f_plus = feval(f,exp(par_plus)); %Solution of ODE at the increment
        f_minus = feval(f,exp(par_minus)); %Solution of ODE at the increment
            
        % Might need to reorient if necessary 
        for j=1:outs
            sens(j,:,i) = (f_plus(j,:) - f_minus(j,:))./(2.*STEP);
        end
    end
end

% If set to 1, normalize the sensitivities with respect to the magnitude of
% the output (e.g., divide by the solution from the non-perturbed
% parameters
if NORMALIZE == 1
    if outs>1
        for k=1:outs
            sens(k,:,:) = sens(k,:,:)./nom_sol(k,:);
        end
    else
        sens = sens./nom_sol';
    end
end
end