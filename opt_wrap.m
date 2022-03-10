% Wrapper file that returns the cost functional of interest

function J = opt_wrap(q,IDS,pars,IC,tspace,NC,data,rflag)
pars(IDS) = q;
J = call_model_DRAM(pars,IC,tspace,NC,data,rflag);
end