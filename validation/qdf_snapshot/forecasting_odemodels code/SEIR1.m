% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function dx=SEIR1(t,x,params0,extra0)

%  params0(1) = beta, params0(2)=k,  params0(3)=gamma, params0(4)=N

dx=zeros(5,1);  % define the vector of the state derivatives

dx(1,1)= -params0(1)*x(1,1).*x(3,1)./params0(4); %S
dx(2,1)= params0(1)*x(1,1).*x(3,1)/params0(4) -params0(2)*x(2,1); %E
dx(3,1)= params0(2)*x(2,1) - params0(3)*x(3,1); %I
dx(4,1)= params0(3)*x(3,1); %R
dx(5,1)= params0(2)*x(2,1); %cumulative infections

