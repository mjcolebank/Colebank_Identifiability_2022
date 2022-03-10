% Routine for obtaining outputs from the model after solving the system of
% ODEs and DAEs

function outputs = get_model_results(t,y,pars)
% Volumes
Vsa = y(1,:);
Vsv = y(2,:);
Vra = y(3,:);
Vrv = y(4,:);
Vpa = y(5,:);
Vpv = y(6,:);
Vla = y(7,:);
Vlv = y(8,:);
Vsw = y(9,:);
% Sarcomere model states
cla = y(10,:);
clv = y(11,:);
cra = y(12,:);
crv = y(13,:);
csw = y(14,:);

Lla = y(15,:);
Llv = y(16,:);
Lra = y(17,:);
Lrv = y(18,:);
Lsw = y(19,:);

% y free wall
ym = y(20,:);



%% Load in parameters
% Parameters for heart chambers
Vpars  = pars.V;
Vwall  = Vpars(1:5);
Am_ref = Vpars(6:10);

V_LVwall = Vwall(2);
V_RVwall = Vwall(4);
V_SWwall = Vwall(5);

% Parameters for the CV models
CVpars = pars.CV;
%% First, solve equations 9-11 to obtain xm, Am, and Cm
% Note, in original Lumen's code we assume that ysw = ylv = yrv; We do a
% second iteration and update eveything else

Vm_LV = - Vlv - 0.5.*(V_LVwall+V_SWwall)+Vsw;
Vm_RV =   Vrv + 0.5.*(V_RVwall+V_SWwall)+Vsw;

Vm         = [Vm_LV;Vm_RV;Vsw];
c          = [clv;crv;csw];
L          = [Llv;Lrv;Lsw];
Vwall_ven  = Vwall([2 4 5]);
Am_ref_ven = Am_ref([2 4 5]);
%% Second iteration: recalculate
%% Now that we have the values of Am and Cm, we want to calculate the
%  mechanics of the curved cardiac walls
[xm,Am,Cm] = get_wallsegment(Vm,ym);

[Tm,sig_f,c,L,eps_f] = solve_wall(t,Vwall_ven,Am,Cm,Am_ref_ven,c,L,pars,0);
[Tx,Ty] = get_subtensions(Tm,xm,ym);


%%
p_trans_LV = 2.*(Tx(1,:)./ym);
p_trans_RV = 2.*(Tx(2,:)./ym);

plv = -p_trans_LV;
prv =  p_trans_RV;


%% Now calculate variables in the atria

% Calculate midwall volume (see Lumen's code)
Vm_la = 0.5*Vwall(1)+Vla;
Vm_ra = 0.5*Vwall(3)+Vra;

% Calulate atrial curvature
Cm_la = (4.*pi./3./Vm_la).^(1/3);
Cm_ra = (4.*pi./3./Vm_ra).^(1/3);

% Now determine the area of the midwall
Am_la = 4.*pi./(Cm_la.^2);
Am_ra = 4.*pi./(Cm_ra.^2);

% Now get the tensions 
atr_offset = pars.SarcA(end);
[Tm_la,G_tot_la,cla,Lla,eps_f_la] = solve_wall(t+atr_offset,Vwall(1),Am_la,Cm_la,Am_ref(1),cla,Lla,pars,1);
[Tm_ra,G_tot_ra,cra,Lra,eps_f_ra] = solve_wall(t+atr_offset,Vwall(3),Am_ra,Cm_ra,Am_ref(3),cra,Lra,pars,1);

% Finally, determine the transmural pressure (different than ventricles,
% see Lumens 2009 code)
pla = abs(2.*Tm_la.*Cm_la);
pra = abs(2.*Tm_ra.*Cm_ra);

%% Now, consider the rest of the circulatory model
Ra_val = CVpars(1);
Rm_val = CVpars(2);
Rp_val = CVpars(3);
Rt_val = CVpars(4);

Rvc    = CVpars(5);
Rpv    = CVpars(6);

Rs     = CVpars(7);
Rp     = CVpars(8);

Csa    = CVpars(9);
Csv    = CVpars(10);
Cpa    = CVpars(11);
Cpv    = CVpars(12);

psa = Vsa/Csa; % systemic arteries 
psv = Vsv/Csv; % systemic veins 
ppa = Vpa/Cpa; % pressure pulmonary artery
ppv = Vpv/Cpv; % pressure pulmonary vein

% right heart valves
qt_val = max((pra-prv)./Rt_val,0); % tricuspid valve opens 
qp_val = max((prv-ppa)./Rp_val,0); % pulmonic valve opens

% left heart valves
qm_val = max((pla-plv)./Rm_val,0); % mitral valve opens
qa_val = max((plv-psa)./Ra_val,0); % aortic valve opens

% venous flow
qvc = max((psv-pra)./Rvc,0);
qpv = (ppv-pla)./Rpv;

% Non-valve systemic/pulmonary flows
qs   = (psa - psv)/Rs;
qp   = (ppa - ppv)/Rp;

% Save to output
outputs.t      = t;
outputs.V      = y(1:9,:);
outputs.Sarc   = y(10:19,:);
outputs.p      = [psa;psv;pra;prv;ppa;ppv;pla;plv];
outputs.q      = [qa_val;qs;qvc;qt_val;qp_val;qp;qpv;qm_val];
outputs.Am     = [Am_la; Am(1,:); Am_ra; Am(2:3,:)];
outputs.Cm     = [Cm_la; Cm(1,:); Cm_ra; Cm(2:3,:)];
outputs.ef     = [eps_f_la; eps_f(1,:); eps_f_ra; eps_f(2:3,:)];
outputs.stress = [G_tot_la; sig_f(1,:); G_tot_ra; sig_f(2:3,:)];


end