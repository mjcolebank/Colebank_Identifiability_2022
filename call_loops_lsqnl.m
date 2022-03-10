%%

function [likelihood,par_set] = call_loops_lsqnl(rflag,pars,par_fix,pars0,par_ids,upper,lower,init,n,data,IC,sc,tspan,id,options)
likelihood = zeros(1,n);
par_set    = zeros(length(init),n);
init_i = init(par_fix);
upper_fix =  upper(par_fix);
lower_fix =  lower(par_fix);
    for j=1:(n/2)
        j
        q_fix = [id,pars(n/2 + j)];
        [x,rsnorm,fval] = lsqnonlin(@(q)call_model_lsqnon(q, q_fix,pars0,par_ids,data,IC,sc,tspan,rflag), init_i,lower_fix, upper_fix, options);
        likelihood(n/2 + j) = sum(fval.^2);
        par_set([id, par_fix],n/2 + j) = [q_fix(2) x];
        init_i = x;
    end
init_i = init(par_fix);
    for j=1:(n/2)
        j
        q_fix = [id,pars(n/2 - (j-1))];
        [x,rsnorm,fval] = lsqnonlin(@(q)call_model_lsqnon(q, q_fix,pars0,par_ids,data,IC,sc,tspan,rflag), init_i,lower_fix, upper_fix, options);
        likelihood(n/2 - (j-1)) = sum(fval.^2);
        par_set([id par_fix],n/2 - (j-1)) = [q_fix(2) x];
        init_i = x;
    end

end