% Function that takes in the the sarcomere fiber strain and computes the
% total stress as a function of both active and pass stresses.

function [G_f_total,dC_dt,dLsc_dt] = sarcomere(t,eps_f,Lsc,C,pars,AV)
if AV==1
    Spars = pars.SarcA;
else
    Spars = pars.SarcV;
end
T     = pars.T;
tmod = mod(t,T);


Ls_ref  = Spars(1);
Ls_iso  = Spars(2);
Vmax    = Spars(3);
Lsc0    = Spars(4);
C_rest  = Spars(5);

% Timing parameters are treated as parameters here: published code uses
% timing relative to some activation time
tauR  = Spars(6);
tauD  = Spars(7);
tauSC = Spars(8);

% Active and passive fiber stress values
sig_act = Spars(9);
sig_pas = Spars(10);

% Parameter describing the roles of titin and ECM
Ls_ref_pas   = Spars(11);
Ls_pas_stiff = Spars(12);
k1           = Spars(13);



% First, calculate the sarcomere length as a function of Sarcomere strain
Ls = Ls_ref.*exp(eps_f); %eq. B1

% Next, define the differential equation for the contractile element
% length, Lsc
dLsc_dt = ((Ls-Lsc)./Ls_iso - 1).*Vmax; %eq. B2

%% Active force
% Now, define the terms needed for the differential equation for sarcomere
% mechanical activiation, C

CL = tanh(4.0.*(Lsc-Lsc0).^2); %eq. B4

xF = min(8,max(0,tmod/tauR)); % eq. B5b
F_rise = 0.02.*xF.^3 .*(8-xF).^2 .* exp(-xF); %eq. B5

T_Lsc = tauSC.*(0.29+0.3.*(Lsc)); %eq. B6

% NOTE: This function has a discontinuity at t=T. Lumens actually uses an
% approximation (not included in the text;
exp_func = 1.0+exp((T_Lsc - tmod)./tauD); %in eq. B3
dC_dt = CL.*F_rise./tauR + ((C_rest-C)./exp_func)./tauD; %eq. B3
%

% Calculate the active stress
G_act  = sig_act.*C.*(Lsc-Lsc0).*(Ls-Lsc)./Ls_iso;

%% Passive force (titin and ECM)
% Now determine the passive stress, which depends on the contributions of
% titin and collagen to passive myocardial stiffness
% this component is based on the work from Vas Osta et al. 2020 and
% Walmsley et al. 2015

%Stretch in the passive components of the wall
Lambda_pas = Ls./Ls_ref_pas;

% First, determine the stress due to ECM
G_ecm = 0.0349.*sig_pas.*(Lambda_pas.^k1 - 1.0);

% Next, determine the stress due to titin, which also depends on active
% force
k2 = 2.0.*(Ls_ref./Ls_pas_stiff);
G_titin = 0.01.*sig_act.*(Lambda_pas.^k2 - 1.0);

G_f_total = G_ecm+G_titin+G_act;

end