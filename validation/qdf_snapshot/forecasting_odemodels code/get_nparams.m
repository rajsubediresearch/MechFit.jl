function numparams=get_nparams(method1,dist1,params_num,fixI0)

numparams=params_num;

numparams=numparams+fixI0;


if method1==0 % Least squares: the Gaussian likelihood behind the AICc formula below
               % estimates the variance (SSE/n), so it counts as one parameter whatever
               % dist1 the bootstrap uses (0, 1 or 2). dist1 only sets the bootstrap noise.

    numparams=numparams+1;

elseif method1==3 | method1==4 %Neg. Binomial requires one more parameter (alpha)

    numparams=numparams+1;

elseif method1==5

    numparams=numparams+2;  %Neg. Binomial requires 2 more parameters (alpha,d)

elseif method1==6 % Laplace distribution for SAD method -- one parameter for the scale parameter

    numparams=numparams+1;

end


