% Function that provides the nominal parameter values for the
% Triseg/circulation model.
% Author: MJ Colebank
% Last edited: 3/10/2022
%
% THIS CODE COMES WITH NO GUARANTEES
% Abbreviations: LA -> left atrium, LV -> left ventricle, RA -> right
% atrium, RV -> right ventricle, S -> septum, Am -> midwall area, Cm ->
% midwall curvature

function [IC,pars] = get_pars_newsarc(data) %make variable input in future versions
% Set the nominal parameter values
T = 60./data.HR;    % Cardiac cycle length (s)

% Riches: 84.7 mL/Kg
Vtot = data.BW*84.7/1000; %microliters

%% Define conversion factors
p_conv = 0.133322;% mmHg -> Kpa
CO_ml = data.CO/60; % Change to ml/s
%% First, CV model parameters
% Relative differences in pressure are based on Boron and Boulapep, Tables 19-3 (vascular) and 22-3 (heart)

% Baseline pressures for parameter estimates
p_lv_sys  = data.p_LV_sys*p_conv;
p_lv_dias = data.p_LV_dias*p_conv;

p_rv_sys  = data.p_RV_sys*p_conv;
p_rv_dias = data.p_RV_dias*p_conv;

p_sa_sys  = 0.99*p_lv_sys;
p_sa_dias = 0.66*p_lv_sys;
p_sa_mean = (p_sa_sys+2.*p_sa_dias)./3;

p_pa_sys  = 0.99*p_rv_sys;
p_pa_dias = 0.27*p_rv_sys;
p_pa_mean = (p_pa_sys+2.*p_pa_dias)./3;

p_sys_cap = 0.26*p_sa_mean;
p_pul_cap = 0.66*p_pa_mean;

p_sv_mean = 0.16*p_sa_mean;
p_pv_mean = 0.33*p_pa_mean;

p_la_sys  = 1.5*p_lv_dias;
p_la_dias = 0.9*p_pv_mean;

p_ra_sys  = 1.25*p_rv_dias;
p_ra_dias = 0.9*p_sv_mean;


%% Now, initialize the volume estimates
% These are the UNS
% Using Lumens 2009 and Boron pg. 878
Vsa_tot   = 0.14*Vtot;   % mL
Vsv_tot   = 0.70*Vtot;   % mL
Vra_tot   = 0.010*Vtot;  % mL
Vrv_tot   = 0.026*Vtot;  % mL
Vpa_tot   = 0.026*Vtot;  % mL
Vpv_tot   = 0.062*Vtot;  % mL
Vla_tot   = 0.010*Vtot;  % mL
Vlv_tot   = 0.026*Vtot;  % mL
Vsw_tot   = 0.0084*Vtot; % mL

Vm_SW = Vsw_tot;  % mL
ym    = 0.262957; % Initial guess for midwall junction (updated later)

%% Calculate stressed volumes based on Beneken and DeWitt 1966

Vsa = Vsa_tot.*0.27;
Vsv = Vsv_tot.*0.075;
Vra = Vra_tot.*1.0;
Vrv = Vrv_tot.*1.0;
Vpa = Vpa_tot.*0.58;
Vpv = Vpv_tot.*0.11;
Vla = Vla_tot.*1.0;
Vlv = Vlv_tot.*1.0;
Vsw = Vsw_tot.*1.0;

%% Next, define TriSeg parameters

% Wall volumes
V_LAwall = 0.0100; % cm^3
V_LVwall = 0.0710; % cm^3
V_RAwall = 0.0100; % cm^3
V_RVwall = 0.0284; % cm^3
V_SWwall = 0.0383; % cm^3

% Reference areas
Am_ref_LA =  0.35; % cm^2
Am_ref_LV =  0.59; % cm^2
Am_ref_RA =  0.35; % cm^2
Am_ref_RV =  0.50; % cm^2
Am_ref_SW =  0.10; % cm^2


%% Parameters for the sarcomere model

% % Atria
Ls_ref  = 2.0;              % micrometer
Ls_iso  = 0.04;             % micrometer
vmax    = 30;               % micrometer per second
Lsc0    = 1.51;             % micrometer
C_rest  = 0.02;             % dimensionless
tauR    = 0.0375*T;         % seconds
tauD    = 0.005*T;          % seconds
tauSC   = 0.15*T;           % seconds
sig_act = 1.5*data.sig_act; % KPa
sig_pas = 1.5*data.sig_pas; % KPa
Ls_ref_pas   = 1.8;         % micrometer
Ls_pas_stiff = 0.6;         % micrometer
k1           = 10;          % dimensionless
t_atr_offset = 0.18.*T;     % seconds

Spars_A = [Ls_ref; Ls_iso; vmax; Lsc0; C_rest; ...
         tauR; tauD; tauSC; sig_act; sig_pas; ...
         Ls_ref_pas; Ls_pas_stiff; k1; t_atr_offset];

% Ventricles
Ls_ref  = 2.0;                % micrometer
Ls_iso  = 0.04;               % micrometer
vmax    = 12.36;              % micrometer per second
Lsc0    = 1.51;               % micrometer
C_rest  = 0.02;               % dimensionless
tauR    = 0.02*T;             % seconds
tauD    = data.relax_fac/1000;% seconds
tauSC   = 0.45*T;             % seconds
sig_act = data.sig_act;       % KPa
sig_pas = data.sig_pas;       % KPa
Ls_ref_pas   = 1.8;           % micrometer
Ls_pas_stiff = 0.6;           % micrometer
k1           = 10;            % dimensionless

Spars_V = [Ls_ref; Ls_iso; vmax; Lsc0; C_rest; ...
         tauR; tauD; tauSC; sig_act; sig_pas; ...
         Ls_ref_pas; Ls_pas_stiff; k1];

%% Nominal resistance values
Ra_val = (p_lv_sys - p_sa_sys)./CO_ml;   % KPa s / muL, Aortic Valve Resistance
Rm_val = (p_la_sys - p_lv_dias)./CO_ml;  % KPa s / muL, Mitral Valve Resistance
Rp_val = (p_rv_sys - p_pa_sys)./CO_ml;   % KPa s / muL, Pulmonic Valve Resistance
Rt_val = (p_ra_sys - p_rv_dias)./CO_ml;  % KPa s / muL, Tricuspid Valve Resistance
Rvc    = (p_sv_mean - p_ra_dias)./CO_ml; % KPa s / muL, Vena Cava Resistance
Rpv    = (p_pv_mean - p_la_dias)./CO_ml; % KPa s / muL, Pulmonary venous Resistance
Rs     = (p_sa_mean-p_sys_cap)./CO_ml;   % KPa s / muL, Systemic vascular Resistance
Rp     = (p_pa_mean-p_pul_cap)./CO_ml;   % KPa s / muL, Pulmonary vascular Resistance
Csa    = Vsa./p_sa_sys;                  % muL / KPa, Systemic artery Compliance
Csv    = Vsv./p_sv_mean;                 % muL / KPa, Systemic venous Compliance
Cpa    = Vpa./p_pa_sys;                  % muL / KPa, Pulmonary artery Compliance
Cpv    = Vpv./p_pv_mean;                 % muL / KPa, Pulmonary venous Compliance


% Set parameter values
Vpars = [V_LAwall; V_LVwall; V_RAwall; V_RVwall; V_SWwall; ...
         Am_ref_LA; Am_ref_LV; Am_ref_RA; Am_ref_RV; Am_ref_SW];


CVpars = [Ra_val; Rm_val; Rp_val; Rt_val; Rvc; Rpv; Rs; Rp; ...
          Csa; Csv; Cpa; Cpv];
      
if any(CVpars<0)
    error('Negative Parameters');
end
      
pars.V     = Vpars;
pars.SarcA = Spars_A;
pars.SarcV = Spars_V;
pars.CV    = CVpars;
pars.T     = T;

%% 
% To get initial values for the septal location states and the sarcomere
% states, we need to solve for zero tension to begin with.
Vm_LV = - Vlv - 0.5.*(V_LVwall+V_SWwall)+Vsw;
Vm_RV =   Vrv + 0.5.*(V_RVwall+V_SWwall)+Vsw;
Vm_LA = 0.5*V_LAwall+Vla;
Vm_RA = 0.5*V_RAwall+Vra;
Vm = [Vm_LA Vm_LV Vm_RA Vm_RV Vm_SW];
[xm,Am,Cm] = get_wallsegment(Vm,ym);

cla0 = C_rest;
clv0 = C_rest;
cra0 = C_rest;
crv0 = C_rest;
csw0 = C_rest;

Lla0 = Ls_ref;
Llv0 = Ls_ref;
Lra0 = Ls_ref;
Lrv0 = Ls_ref;
Lsw0 = Ls_ref;

ym0 = ym;

% Initial guess of initial conditions
IC = [Vsa; Vsv; Vra; Vrv; Vpa; Vpv; Vla; Vlv; Vsw;...
      cla0; clv0; cra0; crv0; csw0; ...
      Lla0; Llv0; Lra0; Lrv0; Lsw0; ...
      ym0];

% Update the septal midwall volume and junction point by solving a root
% finding problem
Vm = Vm([2 4 5]);
[ym0,Vm,~,~] = get_initial_conditions(0,xm,Vm,ym,IC,pars,[Vlv;Vrv]);
IC(9)  = Vm(3);
IC(20) = ym0;

end
