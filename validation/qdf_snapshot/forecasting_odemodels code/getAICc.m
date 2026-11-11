function [AICc,part1,part2,numparams]=getAICc(method1,dist1,params_num,fixI0,fval,n)
  
numparams=params_num;
    
numparams=numparams+fixI0;
    

if method1==0 % Least squares: the Gaussian likelihood behind the AICc formula below
               % estimates the variance (SSE/n), so it counts as one parameter whatever
               % dist1 the bootstrap uses (0, 1 or 2). dist1 only sets the bootstrap noise.

    numparams=numparams+1;

elseif method1==3 | method1==4 %Neg. Binomial requires one more parameter (alpha)

    numparams=numparams+1;

elseif method1==5

    numparams=numparams+2;  % Neg. Binomial requires 2 more parameters (alpha,d)

elseif method1==6  % Laplace distribution -- one parameter for scale

    numparams=numparams+1;

end

% The small-sample correction 2k(k+1)/(n-k-1) is only defined for n - k - 1 > 0,
% i.e. at most n - 2 estimated quantities. Fail clearly instead of returning an
% infinite or negative AICc.
if n-numparams-1 <= 0
    error('getAICc:tooFewData', ...
        'AICc needs at least k + 2 data points: n = %d, k = %d estimated quantities (k must be at most n - 2).', n, numparams);
end

switch method1
    
    case 0
        
        AICc= n*log(fval/n) + n*log(2*pi) + n + 2*numparams + (2*numparams*(numparams+1))/(n-numparams-1);
        
        part1= n*log(fval/n);
        
        part2= n*log(2*pi) + n + 2*numparams + (2*numparams*(numparams+1))/(n-numparams-1);
    
    case 6

        AICc= 2*n*log(2*fval/n) + 2*n + 2*numparams + (2*numparams*(numparams+1))/(n-numparams-1);
        
        part1 = 2*n*log(2*fval/n);
        
        part2=  2*n + 2*numparams + (2*numparams*(numparams+1))/(n-numparams-1);


    otherwise
        
        AICc=-2*(-fval) + 2*numparams + (2*numparams*(numparams+1))/(n-numparams-1);
        
        part1=-2*(-fval);
        
        part2=2*numparams + (2*numparams*(numparams+1))/(n-numparams-1);

end

%method1

%fval

%n

%numparams
