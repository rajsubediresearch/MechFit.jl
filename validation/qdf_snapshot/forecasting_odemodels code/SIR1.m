% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function dx=SIR1(t,x,params0,extra0) %#ok<INUSD>
% SIR model. params0(1) = beta, params0(2) = gamma, params0(3) = N.
% States: S, I, R, C (cumulative infections).
% FIX: the solver always passes four inputs (t, x, params, extra0); the
% shipped version declared only three and failed with "Too many input arguments".

dx=zeros(4,1);  % define the vector of the state derivatives

dx(1,1)= -params0(1)*x(1,1).*x(2,1)./params0(3); %S
dx(2,1)= params0(1)*x(1,1).*x(2,1)./params0(3) - params0(2)*x(2,1); %I
dx(3,1)= params0(2)*x(2,1); %R
dx(4,1)= params0(1)*x(1,1).*x(2,1)./params0(3); %cumulative infections
