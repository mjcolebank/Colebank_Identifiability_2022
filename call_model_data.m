% Wrapper file to call the model
function [y_out,J] = call_model_data(pars,IC,tspace,NC,data)
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

options=odeset('Mass',M,'RelTol',1e-6, 'AbsTol',1e-6);  %sets how accurate the ODE solver is
ysol = ode15s(@DE_model_mass,tspace,IC,options,pars); 

if ysol.x(end)<tspace(end)
    options=odeset('Mass',M,'RelTol',1e-8, 'AbsTol',1e-8);  %sets how accurate the ODE solver is
    ysol = ode15s(@DE_model_mass,[tspace(1) tspace(end)],IC,options,pars);
    if ysol.x(end)<Tend
        options=odeset('Mass',M,'RelTol',1e-12, 'AbsTol',1e-12);  %sets how accurate the ODE solver is
        ysol = ode15s(@DE_model_mass,[tspace(1) tspace(end)],IC,options,pars);
        if ysol.x(end)<Tend
            options=odeset('Mass',M,'RelTol',1e-16, 'AbsTol',1e-16);  %sets how accurate the ODE solver is
            ysol = ode15s(@DE_model_mass,[tspace(1) tspace(end)],IC,options,pars);
            if ysol.x(end)<Tend
                warning('ODE solver stopped prematurely. Recalculating.');
                %%%%%%%%%%%%%%%
                error('Failed');
                %%%%%%%%%%%%%%%
                NC = round(ysol.x(end)./pars.T)-1;
                Tend = pars.T*NC;
                tspace = 0:pars.dt:Tend;
                ysol = ode15s(@DE_model_mass,[tspace(1) tspace(end)],IC,options,pars);
                tspace = 0:pars.dt:Tend;
                nt = length(tspace);
            end
        end
    end
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

% y_out = [p; V; q; Csarc; Lsarc; Am; Cm; ef; stress];
% J = 0;
% return;
% Look at the sensitivity of the residual vector
% Look at the sensitivity of the residual vector
LV_p_sys  = max(p(8,:));
LV_p_dias = min(p(8,:));
LV_V_sys  = min(V(8,:));
LV_V_dias = max(V(8,:));

p_RV = p(4,:);
V_RV = V(4,:);
p_RV = interp1(linspace(0,1,length(p_RV)),p_RV,linspace(0,1,length(data.p)))';
V_RV = interp1(linspace(0,1,length(V_RV)),V_RV,linspace(0,1,length(data.V)))';

r1 = (LV_p_sys-data.p_LV_sys)./data.p_LV_sys;
r2 = (LV_p_dias-data.p_LV_dias)./data.p_LV_dias;
r3 = (LV_V_sys-data.V_LV_sys)./data.V_LV_sys;
r4 = (LV_V_dias-data.V_LV_dias)./data.V_LV_dias;

r5 = (p_RV-data.p)./data.p./sqrt(length(p_RV));
r6 = (V_RV-data.V)./data.V./sqrt(length(V_RV));
% 
% % r5 = (p_RV-data.p)./data.p;
% % r6 = (V_RV-data.V)./data.V;
% 
r7 = (max(V_RV)-max(data.V))./max(data.V);
r8 = (min(V_RV)-min(data.V))./min(data.V);
% 
r9 = (max(p_RV)-max(data.p))./max(data.p);
r10 = (min(p_RV)-min(data.p));%./min(data.p);
% 
EF = (LV_V_dias-LV_V_sys)./LV_V_dias;
% 
% r_full = [r1;r2;r3;r4;r5;r6]; %all data
% % r_full = [r3;r4;r6]; %just volumes
% % r_full = [r3;r4;r7;r8]; %just volumes
r_full = [r1;r2;r3;r4;r5;r6;r7;r8;r9;r10; EF];
% % r_full = [r1; r2; r3; r4; r7; r8; r9; r10]; %just static values
% % y_out = r_full./10;
y_out = r_full;
% % y_out = [EF LV_V_dias LV_V_sys];
% r_full = y_out;
J = r_full'*r_full;%./10;
% y_out = J;
end


function pars = change_pars(pars)
temp.V     = pars(1:10)';
temp.SarcA = pars(11:24)';
temp.SarcV = pars(25:37)';
temp.CV    = pars(38:49)';
temp.T     = pars(50);
pars       = temp;
end

