% Driver file for running the Triseg model used in Colebank & Chesler 2022
% Author: MJ Colebank
% Last edited: 3/10/2022
%
% THIS CODE COMES WITH NO GUARANTEES
% Abbreviations: LA -> left atrium, LV -> left ventricle, RA -> right
% atrium, RV -> right ventricle, S -> septum, Am -> midwall area, Cm ->
% midwall curvature, ef -> wall strain
clear; clc; close all;
%% First, load in nominal parameters and initial conditons
data = MI_Sham_data_dynamic;        % Data is based on the work by J. Philip et al. 2019
[IC,pars] = get_pars_newsarc(data); % Get the initial conditions and parameter values
%% Next, define the time vector
num_cycles = 30;                 % Number of cycles to run the model for
dt         = 1e-4;               % Time stepping for the output and plotting
Tend       = num_cycles.*pars.T; % The end time point for the ODE/DAE solver
tspace     = 0:dt:Tend;          % The time vector
pars.dt    = dt;                 % Append the time stepping

%% Now solve the model
%% MASS MATRIX DAE APPROACH
M = eye(20);  % Define the mass matrix
M(9,9)   = 0; % DAE - Vsw
M(20,20) = 0; % DAE - ym

options=odeset('Mass',M,'RelTol',1e-8, 'AbsTol',1e-8);  %sets how accurate the ODE solver is

%% Timing statements and model solve using ODE15s
tic 
ysol = ode15s(@DE_model_mass,tspace,IC,options,pars); 
toc

%% Recast model solution at specified time points
tic
ysols = deval(ysol,tspace);
toc

%% Obtain all relevant outputs (volumes, pressures, sarcomere dynamics)
n_cyc_plot = 2; %Number of cycles you want to plot
tstart = find(tspace>pars.T*(num_cycles-n_cyc_plot),1);
toutput = tspace(tstart:end);

% Get model outputs
tic
outputs = get_model_results(toutput,ysols(:,tstart:end),pars);
toc


%% Now plot outputs
p      = outputs.p./0.133322;  % KPa  --> mmHg  (see get_model_results for output order)
V      = outputs.V.*1e3;       % mL   --> muL   (see get_model_results for output order)
q      = outputs.q.*1e3;       % mL/s --> muL/s (see get_model_results for output order)
Csarc  = outputs.Sarc(1:5,:);  % Contractility (Gamma in the manuscript) (LA,LV,RA,RV,S)
Lsarc  = outputs.Sarc(6:10,:); % Sarcomere contractile length (Lsc)      (LA,LV,RA,RV,S)
Am     = outputs.Am;           % Midwall area of cardiac chambers        (LA,LV,RA,RV,S)
Cm     = outputs.Cm;           % Curvature of cardiac chambers           (LA,LV,RA,RV,S)
ef     = outputs.ef;           % Strain of the cardiac chamber           (LA,LV,RA,RV,S)
stress = outputs.stress;       % Wall stress of the cardiac chamber      (LA,LV,RA,RV,S)

% Time vector for plotting purposes
t      = tspace(tstart:end) - tspace(tstart);

%% Plotting
figure(1);
subplot(3,3,1); hold on;
plot(t,p(1,:),'LineWidth',3);
ylabel('SA pressure (mmHg)'); set(gca,'FontSize',16);
subplot(3,3,2); hold on;
plot(t,p(2,:),'LineWidth',3);
ylabel('SV pressure (mmHg)'); set(gca,'FontSize',16);
subplot(3,3,3); hold on;
plot(t,p(3,:),'LineWidth',3);
ylabel('RA pressure (mmHg)'); set(gca,'FontSize',16);
subplot(3,3,4); hold on;
plot(t,p(4,:),'LineWidth',3);
ylabel('RV pressure (mmHg)'); set(gca,'FontSize',16);
subplot(3,3,5); hold on;
plot(t,p(5,:),'LineWidth',3);
ylabel('PA pressure (mmHg)'); set(gca,'FontSize',16);
subplot(3,3,6); hold on;
plot(t,p(6,:),'LineWidth',3);
ylabel('PV pressure (mmHg)'); set(gca,'FontSize',16);
subplot(3,3,7); hold on;
plot(t,p(7,:),'LineWidth',3);
ylabel('LA pressure (mmHg)'); set(gca,'FontSize',16);
subplot(3,3,8); hold on;
plot(t,p(8,:),'LineWidth',3);
ylabel('LV pressure (mmHg)'); set(gca,'FontSize',16);


figure(2);
subplot(3,3,1); hold on;
plot(t,V(1,:),'LineWidth',3);
ylabel('SA volume (\mu L)');  set(gca,'FontSize',16);
subplot(3,3,2); hold on;
plot(t,V(2,:),'LineWidth',3);
ylabel('SV volume (\mu L)');  set(gca,'FontSize',16);
subplot(3,3,3); hold on;
plot(t,V(3,:),'LineWidth',3);
ylabel('RA volume (\mu L)'); set(gca,'FontSize',16);
subplot(3,3,4); hold on;
plot(t,V(4,:),'LineWidth',3);
ylabel('RV volume (\mu L)'); set(gca,'FontSize',16);
subplot(3,3,5); hold on;
plot(t,V(5,:),'LineWidth',3);
ylabel('PA volume (\mu L)'); set(gca,'FontSize',16);
subplot(3,3,6); hold on;
plot(t,V(6,:),'LineWidth',3);
ylabel('PV volume (\mu L)'); set(gca,'FontSize',16);
subplot(3,3,7); hold on;
plot(t,V(7,:),'LineWidth',3);
ylabel('LA volume (ml)'); set(gca,'FontSize',16);
subplot(3,3,8); hold on;
plot(t,V(8,:),'LineWidth',3);
ylabel('LV volume (\mu L)'); set(gca,'FontSize',16);
subplot(3,3,9); hold on;
plot(t,V(9,:),'LineWidth',3);
ylabel('SW volume (\mu L)'); set(gca,'FontSize',16);


figure(3);
subplot(2,2,1); hold on;
plot(V(7,:),p(7,:),'LineWidth',3);
title('LA');
set(gca,'FontSize',16);
ylabel('Pressure (mmHg');
subplot(2,2,2); hold on;
plot(V(3,:),p(3,:),'LineWidth',3);
title('RA');
set(gca,'FontSize',16);
subplot(2,2,3); hold on;
plot(V(8,:),p(8,:),'LineWidth',3);
ylabel('Pressure (mmHg');
xlabel('Volume (\mu L)');
title('LV');
set(gca,'FontSize',16);
subplot(2,2,4); hold on;
plot(V(4,:),p(4,:),'LineWidth',3);
title('RV');
xlabel('Volume (\mu L)');
set(gca,'FontSize',16);


figure(4); %clf;
subplot(3,3,1); hold on; plot(ef(2,:),stress(2,:),'LineWidth',3); %xlim([-0.2 0.1]);
xlabel('Myofiber strain'); ylabel('Myofiber stress (KPa)');
set(gca,'FontSize',16); title('LV');grid on;

subplot(3,3,2); hold on;
plot(ef(5,:),stress(5,:),'LineWidth',3); %xlim([-0.2 0.1]);
xlabel('Myofiber strain'); ylabel('Myofiber stress (KPa)');
set(gca,'FontSize',16); title('S');grid on;

subplot(3,3,3); hold on;
plot(ef(4,:),stress(4,:),'LineWidth',3); %xlim([-0.2 0.1]);
xlabel('Myofiber strain'); ylabel('Myofiber stress (KPa)');
set(gca,'FontSize',16); title('RV');grid on;

subplot(3,3,4); hold on;
plot(t,-Cm(2,:),'LineWidth',3); 
xlabel('Time (s)'); ylabel('Midwall curvature (cm^{-1})');
set(gca,'FontSize',16); title('LV');grid on;

subplot(3,3,5); hold on;
plot(t,Cm(5,:),'LineWidth',3); 
xlabel('Time (s)'); ylabel('Midwall curvature (cm^{-1})');
set(gca,'FontSize',16); title('S'); grid on;

subplot(3,3,6); hold on;
plot(t,Cm(4,:),'LineWidth',3);
xlabel('Time (s)'); ylabel('Midwall curvature (cm^{-1})');
set(gca,'FontSize',16); title('RV');grid on;

subplot(3,3,7); hold on;
plot(t,(Am(2,:)-Am(2,1))./Am(2,1),'LineWidth',3); 
xlabel('Time (s)'); ylabel('Engineering Strain');
set(gca,'FontSize',16); title('LV');grid on;

subplot(3,3,8); hold on;
plot(t,(Am(5,:)-Am(5,1))./Am(5,1),'LineWidth',3); 
xlabel('Time (s)'); ylabel('Engineering Strain');
set(gca,'FontSize',16); title('S'); grid on;

subplot(3,3,9); hold on;
plot(t,(Am(4,:)-Am(4,1))./Am(4,1),'LineWidth',3); 
xlabel('Time (s)'); ylabel('Engineering Strain');
set(gca,'FontSize',16); title('RV');grid on;
