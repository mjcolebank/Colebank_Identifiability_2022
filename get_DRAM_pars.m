function params = get_DRAM_pars_NEW(q,IDS,UB,LB)

par_names = {'VW_{la}','VW_{lv}','VW_{ra}','VW_{rv}','VW_{sw}',...
          'Amref_{la}','Amref_{lv}','Amref_{ra}','Amref_{rv}','Amref_{sw}',...
          'Lsref_{a}','Lsiso_{a}','vmax_{a}','Lsc0_{a}','Crest_{a}',...
          'tauR_{a}','tauD_{a}','tauSC_{a}','sigact_{a}','sigpas_{a}',...
          'Lsref-pas_{a}','Ls-pas-stiff{a}','k1_{a}','toffset',...
         'Lsref_{v}','Lsiso_{v}','vmax_{v}','Lsc0_{v}','Crest_{v}',...
          'tauR_{v}','tauD_{v}','tauSC_{v}','sigact_{v}','sigpas_{v}',...
          'Lsref-pas_{v}','Ls-pas-stiff{v}','k1_{v}', ... 
          'R_{a,val}','R_{m,val}','R_{p,val}','R_{t,val}','R_{vc}','R_{pv}',...
          'R_s','R_p','C_{sa}','C_{sv}','C_{pa}','C_{pv}'};


pars = q;

% Construct the "params" structure for DRAM
params = cell(length(IDS),1);
for i=1:length(IDS)
    params{i} = {par_names{IDS(i)},pars(i),LB(i),UB(i)}; %Name
%     params{i,2}   = pars(IDS(i));      %Nominal value
%     params{i,3}   = LB(i);             %Lower bound
%     params{i,4}   = UB(i);             %Lower bound
end


end