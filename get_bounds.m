% Get initial parameter values and upper and lower bound
% Remember to check to see if everything is log-scaled!
function [q0,qall,UB,LB] = get_bounds(q,IDS)
pars = [];
pars(1:10)  = q.V;
pars(11:24) = q.SarcA;
pars(25:37) = q.SarcV;
pars(38:49) = q.CV;
pars(50)    = q.T;

% Now determine the UB and LB
UB = pars(IDS).*5.0;
LB = pars(IDS).*0.05;
q0 = pars(IDS); %

% Change timing parameters to have smaller bounds
% UB(3:7) = 1.2;
% LB(3:7) = 0.8;

q0 = q0.*unifrnd(0.95,1.05,1,length(IDS));
qall = pars;
% qall(IDS) = q0.*pars(IDS);

end
