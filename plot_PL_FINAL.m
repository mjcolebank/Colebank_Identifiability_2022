%% Analyze results from profile likelihood
% This is the final version! 
% MJC 2/28/22
clear; clc; close all;

%%
Names = {'VW_{la}','VW_{lv}','VW_{ra}','VW_{rv}','VW_{s}',...
    'Am_{ref,la}','Am_{ref,lv}','Am_{ref,ra}','Am_{ref,rv}','Am_{ref,s}',...
    'Lsref_{a}','Lsiso_{a}','v_{max,a}','Lsc0_{a}','C_{rest,a}',...
    '\tau_{R,a}','\tau_{D,a}','\tau_{SC,a}','\sigma_{act,a}','\sigma_{pas,a}',...
    'Lsref-pas_{a}','Ls-pas-stiff{a}','k1_{a}','toffset',...
    'Lsref_{v}','Lsiso_{v}','v_{max,v}','Lsc0_{v}','C_{rest,v}',...
    '\tau_{R,v}','\tau_{D,v}','\tau_{SC,v}','\sigma_{act,v}','\sigma_{pas,v}',...
    'Lsref-pas_{v}','Ls-pas-stiff{v}','k1_{v}', ...
    'R_{a,val}','R_{m,val}','R_{p,val}','R_{t,val}','R_{vc}','R_{pv}',...
    'R_s','R_p','C_{sa}','C_{sv}','C_{pa}','C_{pv}'};


plot_colors = [linspace(0.8,0.3,4); linspace(0.8,0.3,4); linspace(0.8,0.3,4)];
for r=1:4
    if r==1
        load likelihood_lsq_R1_Final.mat
    elseif r==2
         load likelihood_lsq_R2_Final.mat
    elseif r==3
        load likelihood_lsq_R3_Final.mat
    else
         load likelihood_lsq_R4_Final.mat
    end
    
    
    
    n = size(par_set,3);
    par_ids = [1:10 13 16:20 22:23 27 30:34 36:37 38:49];
    test_ids = [2 4 5 6 7 9 20 21 22 23 33 34 37];
    par_ids = par_ids(test_ids);
    n_par = length(par_ids);
    n_samp = size(likelihood,2);
    
    
    %% Plot on top
    f_id = 999;
    figure(f_id);
    corr_save = [];
    for i=1:n_par
        q_curr = squeeze(par_set(i,i,:));
        LL = likelihood(i,:);


        ids = find(LL(:)>1e4);
        LL(ids) = [];
        q_curr(ids) = [];

        [LL_min,LL_where] = min(LL);
        LL95 = LL_min+chi2inv(0.95,1);

        id_left = find(LL(1:LL_where)>LL95,1);
        if isempty(id_left)
            x_left = q_curr(1);
        else
            x_left = q_curr(id_left);
        end
        id_right = LL_where + find(LL(LL_where+1:end)>LL95,1);
        if isempty(id_right)
            x_right = q_curr(end);
        else
            x_right = q_curr(id_right);
        end

       
        figure(f_id);
        subplot(4,13,13.*(r-1)+i);hold on;
        plot(q_curr,LL,'-','Color',plot_colors(:,r),'LineWidth',3);
        plot(q_curr(LL_where),LL_min,'ko','MarkerFaceColor','k','MarkerSize',4);
        plot(q_curr,0.*q_curr+(LL_min+chi2inv(0.95,1)),'--','Color',plot_colors(:,r),'LineWidth',3);
        grid on;
        set(gca,'FontSize',10);
%         ylabel('Cost');
%         xlabel(Names{par_ids(i)});
        yticks(LL_min+[chi2inv(0.68,1),chi2inv(0.90,1),chi2inv(0.95,1),chi2inv(0.99,1)]);
        if i==1
        yticklabels({'68%','90%','95%','99%'})
        else
            yticklabels({});
        end
        ylim([LL_min,LL_min+chi2inv(0.991,1)]);
%         xlim([x_left x_right]);
        xlim([0.5 1.5])
        hold off;
        
        
        
        
    end
end
for i=1:13
    subplot(4,13,i); title('');%title(Names(par_ids(i)))
end

