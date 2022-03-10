% Wrapper file to call the model
function r_full = call_model_lsqnon(pars,par_fix,pars0,par_ids,data,IC,sc,tspace,rflag)
if par_fix~=0
    % Only optimize for one parameter at a time
    ID = par_ids(par_fix(1));
    not_IDS = 1:length(sc);
    not_IDS(par_fix(1)) = [];
    par_val  = par_fix(2);
    npar = length(pars);
    pars_full = pars0;
    pars_full(par_ids(not_IDS)) = pars.*sc(not_IDS);
    pars_full(ID) = par_val.*sc(par_fix(1));
    pars = pars_full;
else
    pars = pars0;
end
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
NC = 30;


%%
% NOTE: since we are changing parameters, the IC for the DAE need to be
% recalibrated, especially when changing sarcomere dynamics
%%
Vm_LV = - IC(8) - 0.5.*(pars.V(2)+pars.V(5))+IC(9);
Vm_RV =   IC(4) + 0.5.*(pars.V(4)+pars.V(5))+IC(9);
Vm = [Vm_LV Vm_RV IC(9)];
[xm,Am,Cm] = get_wallsegment(Vm,IC(20));
[ynew,Vm,Tx0,Ty0] = get_initial_conditions(0,xm,Vm,IC(20),IC,pars,[IC(8);IC(4)]);

if abs(Tx0)<1e-6 && abs(Ty0)<1e-6

    IC(9) = Vm(3);
    IC(20) = ynew;
    %% MASS MATRIX DAE APPROACH
    M = eye(20);
    M(9,9)   = 0; %DAE - Vsw
    M(20,20) = 0; %DAE - ym

    nt   = length(tspace);
    Tend = tspace(end);
    pars.dt = 1e-4;

    options=odeset('Mass',M,'RelTol',1e-10, 'AbsTol',1e-10);  %sets how accurate the ODE solver is
    ysol = ode15s(@DE_model_mass,tspace,IC,options,pars);
    Bad_sol = 0;
    if ysol.x(end)<tspace(end)
        warning('ODE solver stopped prematurely. Recalculating.');
        Bad_sol=1;
    end
else
    Bad_sol=1;
end

p_la_data = data.p_la_data(1:8:1083);
p_ra_data = data.p_ra_data(1:8:1083);
p_lv_data = data.p_lv_data(1:8:1083);
p_rv_data = data.p_rv_data(1:8:1083);

V_la_data = data.V_la_data(1:8:1083);
V_ra_data = data.V_ra_data(1:8:1083);
V_lv_data = data.V_lv_data(1:8:1083);
V_rv_data = data.V_rv_data(1:8:1083);

if Bad_sol==1
    if rflag==1
        r_full = p_rv_data(:);
    elseif rflag==2
        r_full = [p_rv_data(:); V_rv_data(:)];
    elseif rflag==3
        r_full = [p_rv_data(:); V_rv_data(:); 0; 0; 0; 0;];
    else
        r_full = [p_rv_data(:); V_rv_data(:); p_lv_data(:); V_lv_data(:)];
    end
    r_full = r_full.*0+1e8;
    return;
end

tstart = max(find(tspace>pars.T*(NC-1),1),find(tspace<=pars.T,1));
toutput = tspace(tstart:nt);
ysols = deval(ysol,toutput);

out = get_model_results(toutput,ysols,pars);
p = out.p./0.133322;%133.322;
V = out.V.*1e3;
temp = V(1:8,:);
if any(p(:)<0) || any(temp(:)<0)
	if rflag==1
		r_full = p_rv_data(:);
	elseif rflag==2
		r_full = [p_rv_data(:); V_rv_data(:)];
	elseif rflag==3
		r_full = [p_rv_data(:); V_rv_data(:); 0; 0; 0; 0];
	else
		r_full = [p_rv_data(:); V_rv_data(:); p_lv_data(:); V_lv_data(:)];
	end
	r_full = r_full.*0+1e8;
	return
end

LV_p_sys  = max(p(8,:));
LV_p_dias = min(p(8,:));
LV_V_sys  = min(V(8,:));
LV_V_dias = max(V(8,:));

RV_p_sys  = max(p(4,:));
RV_p_dias = min(p(4,:));
RV_V_sys  = min(V(4,:));
RV_V_dias = max(V(4,:));

p_LV = p(8,:);
V_LV = V(8,:);
p_RV = p(4,:);
V_RV = V(4,:);

n_data = length(p_lv_data);

p_LV = interp1(linspace(0,1,length(p_LV)),p_LV,linspace(0,1,n_data));
V_LV = interp1(linspace(0,1,length(V_LV)),V_LV,linspace(0,1,n_data));
p_RV = interp1(linspace(0,1,length(p_RV)),p_RV,linspace(0,1,n_data));
V_RV = interp1(linspace(0,1,length(V_RV)),V_RV,linspace(0,1,n_data));

% Dynamic residuals
r_LVP = (p_LV-p_lv_data);
r_LVV = (V_LV-V_lv_data);
r_RVP = (p_RV-p_rv_data);
r_RVV = (V_RV-V_rv_data);

% Static residuals
rS_LVP_sys  = (LV_p_sys-max(p_lv_data));
rS_LVP_dias = (LV_p_dias-min(p_lv_data));

rS_LVV_sys  = (LV_V_sys-min(V_lv_data));
rS_LVV_dias = (LV_V_dias-max(V_lv_data));


if rflag==1
    % Just RV pressure
    r_full = [r_RVP(:)];
elseif rflag==2
    % RV P/V
    r_full = [r_RVP(:); r_RVV(:)];

elseif rflag==3
    % RV dynamic/ LV static
    r_full = [r_RVP(:); r_RVV(:); rS_LVP_sys; rS_LVP_dias; rS_LVV_sys; rS_LVV_dias];
else
    % RV/LV dynamic
    r_full = [r_RVP(:); r_RVV(:); r_LVP(:); r_LVV(:)];
end


end


function pars = change_pars(pars)
temp.V     = pars(1:10)';
temp.SarcA = pars(11:24)';
temp.SarcV = pars(25:37)';
temp.CV    = pars(38:49)';
temp.T     = pars(50);
pars       = temp;
end

