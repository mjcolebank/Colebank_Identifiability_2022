% Analyze the results from the sensitivity analysis
clear; clc; close all
%%
Names = {'VW_{la}','VW_{lv}','VW_{ra}','VW_{rv}','VW_{s}',...
          'Amref_{la}','Amref_{lv}','Amref_{ra}','Amref_{rv}','Amref_{s}',...
          'Lsref_{a}','Lsiso_{a}','vmax_{a}','Lsc0_{a}','Crest_{a}',...
          'tauR_{a}','tauD_{a}','tauSC_{a}','sigact_{a}','sigpas_{a}',...
          'Lsref-pas_{a}','Ls-pas-stiff{a}','k1_{a}','toffset',...
         'Lsref_{v}','Lsiso_{v}','vmax_{v}','Lsc0_{v}','Crest_{v}',...
          'tauR_{v}','tauD_{v}','tauSC_{v}','sigact_{v}','sigpas_{v}',...
          'Lsref-pas_{v}','Ls-pas-stiff{v}','k1_{v}', ... 
          'Ra_val','Rm_val','Rp_val','Rt_val','Rvc','Rpv',...
          'Rs','Rp','Csa','Csv','Cpa','Cpv'};
%%
% load sens_Centered_9_24_2021_out.mat
% load sens_centered_11_18_2021_states.mat
load sens_lognorm_centered_11_19_2021_states.mat


Names = Names(par_ids);


p      = squeeze(sens(1:8,:,:));
V      = squeeze(sens(9:17,:,:));
q      = squeeze(sens(18:25,:,:));
Csarc  = squeeze(sens(26:30,:,:));
Lsarc  = squeeze(sens(31:35,:,:));
Am     = squeeze(sens(36:40,:,:));
Cm     = squeeze(sens(41:45,:,:));
ef     = squeeze(sens(46:50,:,:));
stress = squeeze(sens(51:55,:,:));
% 
[M,~,N] = size(sens); 
sens_norm = zeros(M,N);
% 
for i=1:M
    for j=1:N
        S = sens(i,:,j);
        sens_norm(i,j) = S*S';
    end
end
%% USE FOR RESIDUAL/COST
% load sens_Centered_9_2_2021_data.mat
% load sens_postopt_9_7_2021_EF.mat
% load sens_postopt_9_7_2021_data.mat
% load sens_Centered_9_9_2021_static_parbest.mat
% load sens_Centered_9_13_2021_parbest.mat

% sens = squeeze(sens);
% [M,N] = size(sens); 
% sens_norm = zeros(1,N);
% 
% for j=1:N
%     S = sens(:,j);
%     sens_norm(j) = S'*S;
% end


%%
% outs_ID = 1;
% outs = {'Residual'};
% outs_ID = [4 5 8 11 15 36:40 46:50 51:55];
% outs = {'p_{rv}','p_{pa}','p_{lv}', ...
%         'V_{rv}','V_{lv}', ...
%         'Am_{la}','Am_{lv}','Am_{ra}','Am_{rv}','Am_{sw}',...
%         'ef_{la}','ef_{lv}','ef_{ra}','ef_{rv}','ef_{sw}',...
%         '\sigma_{la}','\sigma_{lv}','\sigma_{ra}','\sigma_{rv}','\sigma_{sw}',...
%         };

outs_ID = [4 8 12 16];
outs = {'p_{rv}','p_{lv}','V_{rv}','V_{lv}'};
Isens_all = [];

% outs_ID = [16];
% outs = {'V_{lv}'};
% Isens_all = [];
% 
% outs_ID = [5];
% outs = {'P_{pa}'};
Isens_all = [];

counter = 1;
for id=outs_ID
% id = 4;
[who,where] = sort(sens_norm(id,:),'descend');
who = who./who(1);
figure(id); clf;
semilogy(who,'o','LineWidth',2); 
xticks(1:N);
xtickangle(45);
xticklabels(Names(where));
ylabel(outs{counter});
set(gca,'FontSize',14);
grid on;
% xlim([1 16]);

% Names(par_ids(where(1:18)))

counter = counter+1;
sens_tol = sqrt(1e-8);%*10;
hold on; plot(1:N,0*who+sens_tol,'--k','LineWidth',3);
Isens_all(end+1,:) = where;
end

%% Try SVD
pars0 = struct2cell(pars);
pars0 = cell2mat(pars0);
pars0 = pars0(par_ids);
for j = outs_ID
    S = squeeze(sens(j,:,:));
    F = S'*S;
    [eU,eS,eV] = svd(F,'econ');
    [Q,R,imap] = qr(eV',0);
    eig_val = diag(eS);
    QR_map(i,:)     = imap;
    [who,where] = find(eig_val./eig_val(1)>(sqrt(1e-8)));
%     dens_store(end+1:end+length(where)) = imap(where);
    figure(1000+j);
    semilogy(diag(eS),'*');
    xticks(1:38);
    xticklabels(Names(imap));
    xtickangle(45)
    set(gca,'FontSize',20);
    Names(imap(who))
    
    V = inv(F);
    SE = sqrt(diag(V));
    [who,where] = sort(SE,'descend');
    figure(2000+j); semilogy(who);
    xticks(1:38);
    xticklabels(Names(where));
    xtickangle(45)
    set(gca,'FontSize',20);
    
    selection_score = norm(sqrt(diag(V))./pars0);
    [who,where] = sort(selection_score,'descend');
    figure(3000+j); semilogy(who,'o');
    xticks(1:38);
    xticklabels(Names(where));
    xtickangle(45)
    set(gca,'FontSize',20);
end
%%
%% Try Eigenvalue method
for j = outs_ID
    S = squeeze(sens(j,:,:));
    F = S'*S;
    [V,D] = eig(F);
    s2vals = diag(D);
    s_temp = min(s2vals);
    par_fix = [];
    par_ids_temp = 1:length(par_ids);
    while s_temp<1e-3
        [who,where] = max(V(:,1));
        par_fix(end+1) = par_ids_temp(where);
        par_ids_temp(where) = [];
        S(:,where) = [];
        F = S'*S;
        [V,D] = eig(F);
        s2vals = diag(D);
        s_temp = min(s2vals);
    end
    
    snorm = sum(S.^2);
    [who,where] = sort(snorm,'descend');
    figure(10000+j); clf;
    semilogy(who,'o','LineWidth',2);
    xticks(1:length(F));
    xtickangle(45);
    xticklabels(Names(par_ids_temp(where)));
    set(gca,'FontSize',14);
    grid on;
    Names(par_fix)
end

