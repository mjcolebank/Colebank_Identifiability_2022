% Wrapper file to call the model
function [y_out,J] = call_model_sens(pars,IC,tspace,NC,data)
if ~isstruct(pars)
    %     pars = exp(pars);
    pars = change_pars(pars);
else
    pars.V     = pars.V;%exp(pars.V);
    pars.SarcA = pars.SarcA;%exp(pars.SarcA);
    pars.SarcV = pars.SarcV;%exp(pars.SarcV);
    pars.CV    = pars.CV;%exp(pars.CV);
    pars.T     = pars.T;%exp(pars.T);
end
%%
% NOTE: since we are changing parameters, the IC for the DAE need to be
% recalibrated, especially when changing sarcomere dynamics
%%
Vm_LV = - IC(8) - 0.5.*(pars.V(2)+pars.V(5))+IC(9);
Vm_RV =   IC(4) + 0.5.*(pars.V(4)+pars.V(5))+IC(9);
Vm = [Vm_LV Vm_RV IC(9)];
[xm,Am,Cm] = get_wallsegment(Vm,IC(20));
[ynew,Vm,~,~] = get_initial_conditions(0,xm,Vm,IC(20),IC,pars,[IC(8);IC(4)]);
IC(9) = Vm(3);
IC(20) = ynew;
%% MASS MATRIX DAE APPROACH
M = eye(20);
M(9,9)   = 0; %DAE - Vsw
M(20,20) = 0; %DAE - ym

nt   = length(tspace);
Tend = tspace(end);
pars.dt = 1e-4;

options=odeset('Mass',M,'RelTol',1e-12, 'AbsTol',1e-12);  %sets how accurate the ODE solver is
ysol = ode15s(@DE_model_mass,tspace,IC,options,pars);

if ysol.x(end)<tspace(end)
    warning('ODE solver stopped prematurely. Recalculating.');
    %%%%%%%%%%%%%%%
    error('Failed');
    %%%%%%%%%%%%%%%

end

ysols = deval(ysol,tspace);

tstart = max(find(tspace>pars.T*(NC-1),1),find(tspace<=pars.T,1));
toutput = tspace(tstart:nt);
out = get_model_results(toutput,ysols(:,tstart:nt),pars);
p = out.p./0.133322;%133.322;
V = out.V.*1e3;
q = out.q.*1e3;
Csarc = out.Sarc(1:5,:);
Lsarc = out.Sarc(6:10,:);
Am = out.Am;
Cm = out.Cm;
ef = out.ef;
stress = out.stress;

y_out = [p; V; q; Csarc; Lsarc; Am; Cm; ef; stress];
J = 0;

end


function pars = change_pars(pars)
temp.V     = pars(1:10)';
temp.SarcA = pars(11:24)';
temp.SarcV = pars(25:37)';
temp.CV    = pars(38:49)';
temp.T     = pars(50);
pars       = temp;
end

