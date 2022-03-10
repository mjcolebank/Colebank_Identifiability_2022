% Wrapper for dram in parallel

function [results,chain,s2chain] = call_dram(q0,F,data,IDS,UB,LB,iter,rflag)

disp([iter q0])
Jopt = @(q) sum(F(q,0));
%% Run an optimization to get the best parameters
options=optimoptions('fmincon', 'Display','iter', ...
     'StepTolerance',  1.0000e-8,'FiniteDifferenceStepSize',1e-4,'FiniteDifferenceType', 'forward',...
     'MaxFunctionEvaluations', 8000, 'MaxIterations', 8000);
[qopt,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(Jopt,q0,[],[],[],[],LB,UB,[],options); 

qopt = q0;


%% DRAM algorithm parameters
method = 'dram';
par_struct = get_DRAM_pars(qopt,IDS,UB,LB);

if rflag==1
    model.N = [136];
    data_DRAM = cell(1,1);
    data_DRAM{1,1}.ydata = data.p_rv_data(1:8:end);
elseif rflag==2
    model.N = [136 136];
    data_DRAM = cell(2,1);
    data_DRAM{1,1}.ydata = data.p_rv_data(1:8:end);
    data_DRAM{2,1}.ydata = data.V_rv_data(1:8:end);
elseif rflag==3
    model.N = [136 136 2 2];
    data_DRAM = cell(4,1);
    data_DRAM{1,1}.ydata = data.p_rv_data(1:8:end);
    data_DRAM{2,1}.ydata = data.V_rv_data(1:8:end);
    data_DRAM{3,1}.ydata = [max(data.p_lv_data) min(data.p_lv_data)];
    data_DRAM{4,1}.ydata = [max(data.V_lv_data) min(data.V_lv_data)];
else
    model.N = [136 136 136 136];
    data_DRAM = cell(4,1);
    data_DRAM{1,1}.ydata = data.p_rv_data(1:8:end);
    data_DRAM{2,1}.ydata = data.V_rv_data(1:8:end);
    data_DRAM{3,1}.ydata = data.p_lv_data(1:8:end);
    data_DRAM{4,1}.ydata = data.V_lv_data(1:8:end);
end
model.ssfun  = F;
J0 = F(qopt,0);
s20 = J0./(sum(model.N)-length(q0));
model.sigma2 = s20; %Sample variance estimate

if cond(HESSIAN)>1e12
    V = eye(length(q0)).*1e-3;
else
    V = mean(s20).*inv(HESSIAN); 
end

DRAM_options.qcov = V; % Sample covariance
DRAM_options.updatesigma = 1;
DRAM_options.method = method;
DRAM_options.nsimu = 50000;
DRAM_options.burnintime = 2000;

[results, chain, s2chain]= mcmcrun(model,data_DRAM,par_struct,DRAM_options);

end
