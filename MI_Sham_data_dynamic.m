% Generate a data file with data necessary for the MI mice
% Author: MJ Colebank
% Last edited: 3/10/2022
%
% THIS CODE COMES WITH NO GUARANTEES
% This data is based on the study by Philip et al. 2019, doi:10.1152/ajpheart.00319.2018
function data = MI_Sham_data_dynamic

data.BW = 29.719;      % g
data.HR = 554;         % in bpm

data.p_LV_dias = 2.0;  %mmHg
data.p_LV_sys  = 76.0; %mmHg

data.p_RV_dias = 0.74; %mmHg
data.p_RV_sys  = 19.6; %mmHg

data.V_LV_sys = 41.40;  %microlieters
data.V_LV_dias = 86.77; %microlieters

data.V_RV_sys = 18.2;   %microliters
data.V_RV_dias = 37.3;  %microliters

data.CO = 10.6;         %in mL per min

data.relax_fac = 8.9;   %milliseconds

% Sarcomere parameters
% These come from Wang et al. 2018, doi: 10.1152/japplphysiol.00725.2017
data.sig_act = 25.52;   %Kpa
data.sig_pas = 1.49;    %Kpa

end